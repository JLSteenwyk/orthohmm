from benchmark_tools.convert_orthomcl_blast import convert_blast, read_fasta_lengths


def test_read_fasta_lengths(tmp_path):
    fasta = tmp_path / "all.fa"
    fasta.write_bytes(b">query description\nAAA\nAA\n>subject\nMMMM\n")

    assert read_fasta_lengths(fasta) == {b"query": 5, b"subject": 4}


def test_convert_blast_matches_orthomcl_bpo_layout(tmp_path):
    fasta = tmp_path / "all.fa"
    blast = tmp_path / "all.blast"
    output = tmp_path / "all.bpo"
    fasta.write_bytes(b">query\n" + b"A" * 20 + b"\n>subject\n" + b"M" * 30 + b"\n")
    blast.write_bytes(
        b"query\tsubject\t80.00\t10\t2\t0\t1\t10\t20\t29\te-20\t100\n"
        b"query\tsubject\t50.00\t8\t4\t0\t11\t18\t3\t8\t2e-10\t50\n"
        b"subject\tsubject\t100.00\t30\t0\t0\t1\t30\t1\t30\t0.0\t200\n"
    )

    count = convert_blast(blast, fasta, output, progress_every=0)

    assert count == 2
    assert output.read_bytes() == (
        b"1;query;20;subject;30;1e-20;59;1:1-10:20-29.2:11-18:3-8.\n"
        b"2;subject;30;subject;30;0.0;100;1:1-30:1-30.\n"
    )
