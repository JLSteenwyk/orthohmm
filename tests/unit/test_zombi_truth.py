import pytest

from benchmark_tools.zombi_truth import event_graph, ortholog_truth, read_fasta, xml_graph, validate_run


def test_speciation_followed_by_duplication_retains_coorthologs():
    graph = {"Root_1": ("S", ("A_1", "B_1")), "A_1": ("D", ("A_2", "A_3")),
             "A_2": ("F", ()), "A_3": ("F", ()), "B_1": ("F", ())}
    assert ortholog_truth(graph) == ({"A_2", "A_3", "B_1"}, {("A_2", "B_1"), ("A_3", "B_1")})


def test_ancestral_duplication_does_not_create_cross_copy_orthologs():
    graph = {"Root_1": ("D", ("Root_2", "Root_3")),
             "Root_2": ("S", ("A_1", "B_1")), "Root_3": ("S", ("A_2", "B_2")),
             **{n: ("F", ()) for n in ("A_1", "A_2", "B_1", "B_2")}}
    assert ortholog_truth(graph)[1] == {("A_1", "B_1"), ("A_2", "B_2")}
    graph["B_2"] = ("L", ())
    assert ortholog_truth(graph) == ({"A_1", "A_2", "B_1"}, {("A_1", "B_1")})


@pytest.mark.parametrize("graph", [
    {"Root_1": ("T", ())}, {"Root_1": ("S", ("missing", "other"))},
    {"Root_1": ("D", ("Root_1", "A_1")), "A_1": ("F", ())},
    {"Root_1": ("F", ()), "extra": ("F", ())},
    {"Root_1": ("S", ("A_1", "A_2")), "A_1": ("F", ()), "A_2": ("F", ())},
])
def test_invalid_histories_fail_closed(graph):
    with pytest.raises(ValueError):
        ortholog_truth(graph)


def test_event_table_validation(tmp_path):
    path = tmp_path / "events.tsv"
    valid = "TIME\tEVENT\tNODES\n0\tO\tRoot\n1\tS\tRoot;1;A;2;B;3\n2\tF\tA;2\n2\tL\tB;3\n"
    path.write_text(valid)
    assert ortholog_truth(event_graph(path)) == ({"A_2"}, set())
    for bad in (valid.replace("\tS\t", "\tT\t"), valid + "3\tF\tA;2\n",
                valid.replace("2\tL", "0.5\tL"), valid.replace("0\tO\tRoot", "0\tO\tA")):
        path.write_text(bad)
        with pytest.raises(ValueError):
            event_graph(path)


def test_xml_event_and_location_validation(tmp_path):
    path = tmp_path / "tree.xml"
    text = '<recGeneTree><phylogeny><clade><name>Root_1</name><eventsRec><loss speciesLocation="Root"/></eventsRec></clade></phylogeny></recGeneTree>'
    path.write_text(text)
    assert xml_graph(path) == {"Root_1": ("L", ())}
    for bad in (text.replace('speciesLocation="Root"', 'speciesLocation="other"'), text.replace("<loss", "<transfer")):
        path.write_text(bad)
        with pytest.raises(ValueError):
            xml_graph(path)


def test_protein_ids_and_alphabet(tmp_path):
    path = tmp_path / "sequences.fasta"
    path.write_text(">a\nACDE\n")
    assert read_fasta(path) == {"a": "ACDE"}
    for text in (">a\nACDE\n>a\nACDE\n", ">a\nAX*\n", ">a\n"):
        path.write_text(text)
        with pytest.raises(ValueError):
            read_fasta(path)


def fixture_run(root):
    files = {
        "T/ExtantTree.nwk": "(A:1,B:1)Root;",
        "G/Genomes/A_GENOME.tsv": "POSITION\tGENE_FAMILY\tORIENTATION\tGENE_ID\n0\t1\t+\t2\n",
        "G/Genomes/B_GENOME.tsv": "POSITION\tGENE_FAMILY\tORIENTATION\tGENE_ID\n0\t1\t+\t3\n",
        "G/Gene_families/1_events.tsv": "TIME\tEVENT\tNODES\n0\tO\tRoot\n1\tS\tRoot;1;A;2;B;3\n2\tF\tA;2\n2\tF\tB;3\n",
        "G/Gene_trees/1_prunedtree.nwk": "(A_2:1,B_3:1)Root_1;",
        "G/Gene_trees/1_rec.xml": '<recGeneTree><phylogeny><clade><name>Root_1</name><eventsRec><speciation speciesLocation="Root"/></eventsRec><clade><name>A_2</name><eventsRec><P speciesLocation="A"/></eventsRec></clade><clade><name>B_3</name><eventsRec><P speciesLocation="B"/></eventsRec></clade></clade></phylogeny></recGeneTree>',
        "S/1_complete.fasta": ">Root_1\nAAAA\n>A_2\nACDE\n>B_3\nFGHI\n",
    }
    for relative, text in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)


def test_complete_crosscheck_excludes_ancestors(tmp_path):
    fixture_run(tmp_path)
    report, sequences = validate_run(tmp_path)
    assert report["extant_genes"] == 2
    assert report["ortholog_pairs"] == [("F1__A_2", "F1__B_3")]
    assert sequences == {"F1__A_2": ("A", "ACDE"), "F1__B_3": ("B", "FGHI")}


@pytest.mark.parametrize("artifact,old,new", [
    ("G/Gene_trees/1_rec.xml", "speciation", "duplication"),
    ("G/Gene_trees/1_prunedtree.nwk", "B_3", "B_9"),
    ("S/1_complete.fasta", ">B_3", ">B_9"),
    ("G/Genomes/B_GENOME.tsv", "\t3\n", "\t9\n"),
])
def test_inconsistent_simulator_outputs_rejected(tmp_path, artifact, old, new):
    fixture_run(tmp_path)
    path = tmp_path / artifact
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError):
        validate_run(tmp_path)
