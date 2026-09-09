import pytest

from benchmark_tools.parallelize_orthomcl_pairs import parallelize_source


SERIAL_SOURCE = """before
for(my $i=0;$i<scalar(@taxa)-1;$i++) {
\tfor(my $j=$i+1;$j<scalar(@taxa);$j++) {
\t\twrite_log("$taxa[$i] and $taxa[$j]");
\t\t$connect{$taxa[$i].' '.$taxa[$j]} = [$taxa[$i], $taxa[$j]];
\t}
}

%blastquery=();
after
"""


def test_parallelizes_original_pair_body():
    result = parallelize_source(SERIAL_SOURCE)

    assert "ORTHOMCL_PAIR_WORKERS" in result
    assert 'write_log("$ta and $tb")' in result
    assert "$process_intertaxon_pair->(@$pair)" in result
    assert "Storable::nstore" in result
    assert "%blastquery=();" in result


def test_rejects_second_transformation():
    with pytest.raises(ValueError, match="already pair-parallel"):
        parallelize_source(parallelize_source(SERIAL_SOURCE))
