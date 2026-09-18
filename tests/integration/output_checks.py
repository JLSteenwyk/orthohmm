"""Order-independent output checks against unchanged biological fixtures."""

from collections import Counter

from Bio import SeqIO


def group_id(label):
    value = label.removesuffix(":")
    assert value.startswith("OG") and value[2:].isdigit()
    return "OG" + str(int(value[2:]))


def groups(path):
    result, seen = {}, set()
    for line in path.read_text().splitlines():
        label, *genes = line.split()
        key = group_id(label)
        assert genes and key not in result
        assert len(set(genes)) == len(genes) and not seen.intersection(genes)
        result[key] = frozenset(genes)
        seen.update(genes)
    return result


def fasta(path):
    result = {}
    for entry in SeqIO.parse(path, "fasta"):
        assert entry.id not in result
        assert entry.description == entry.id
        result[entry.id] = str(entry.seq)
    assert result
    return result


def check_outputs(output, expected, input_files):
    source, species = {}, {}
    for path in input_files:
        for entry in SeqIO.parse(path, "fasta"):
            assert entry.id not in source
            source[entry.id] = str(entry.seq)
            species[entry.id] = path.name
    observed = groups(output / "orthohmm_orthogroups.txt")
    truth = groups(expected / "orthohmm_orthogroups.txt")
    assert set(observed.values()) == set(truth.values())
    assert set().union(*observed.values()) == set(source)
    lines = (output / "orthohmm_gene_count.txt").read_text().splitlines()
    header, *columns = lines[0].split()
    assert header == "files:" and len(set(columns)) == len(columns)
    assert set(columns) == {path.name for path in input_files}
    counts = {}
    for line in lines[1:]:
        label, *values = line.split()
        key = group_id(label)
        assert key not in counts and len(values) == len(columns)
        row = dict(zip(columns, map(int, values)))
        actual = Counter(species[gene] for gene in observed[key])
        assert row == {name: actual[name] for name in columns}
        counts[key] = row
    assert set(counts) == set(observed)
    all_files = list((output / "orthohmm_orthogroups").iterdir())
    assert {path.name for path in all_files} == {key + ".fa" for key in observed}
    for key, genes in observed.items():
        assert fasta(output / "orthohmm_orthogroups" / (key + ".fa")) == {g: source[g] for g in genes}
    # Membership and occupancy independently determine the intended single-copy set.
    selected = {key for key, genes in observed.items()
                if len({species[g] for g in genes}) == len(genes)
                and len(genes) / len(columns) > .5}
    names = (output / "orthohmm_single_copy_orthogroups.txt").read_text().splitlines()
    assert len(names) == len(set(names)) and set(names) == selected
    single_files = list((output / "orthohmm_single_copy_orthogroups").iterdir())
    assert {path.name for path in single_files} == {key + ".fa" for key in selected}
    for key in selected:
        target = {}
        for gene in observed[key]:
            name = species[gene]
            for extension in (".fa", ".faa", ".fas", ".fasta", ".pep", ".prot"):
                if name.endswith(extension):
                    name = name[:-len(extension)]
                    break
            target[name + "|" + gene] = source[gene]
        assert fasta(output / "orthohmm_single_copy_orthogroups" / (key + ".fa")) == target
    return len(observed), len(selected)
