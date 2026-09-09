import pytest

from benchmark_tools.parallelize_orthomcl_pairs import parallelize_source


SERIAL_SOURCE = """before
for(my $i=0;$i<scalar(@taxa)-1;$i++) {
\tfor(my $j=$i+1;$j<scalar(@taxa);$j++) {
\t\tfor(my $k=0;$k<scalar(@nodes1);$k++) {
\t\t\tfor(my $l=0;$l<scalar(@nodes2);$l++) {
\t\t\t\t\tif (blastqueryab($nodes1[$k],$nodes2[$l])) {
\t\t\t\t\t\tmy ($s,$pm,$pe,$pi)=(blastqueryab($nodes1[$k],$nodes2[$l]))[0,3,4,5];
\t\t\t\t\t}
\t\t\t\t\tif (blastqueryab($nodes2[$l],$nodes1[$k])) {
\t\t\t\t\t\tmy ($s,$pm,$pe,$pi)=(blastqueryab($nodes2[$l],$nodes1[$k]))[0,3,4,5];
\t\t\t\t\t}
\t\t\t}
\t\t}
\t\t$connect{$taxa[$i].' '.$taxa[$j]} = [$taxa[$i], $taxa[$j]];
\t}
}

%blastquery=();
after
"""


def test_parallelizes_original_pair_body():
    result = parallelize_source(SERIAL_SOURCE)

    assert "ORTHOMCL_PAIR_WORKERS" in result
    assert "$process_intertaxon_pair->(@$pair)" in result
    assert "Storable::nstore" in result
    assert "Reusing completed inter-taxon result" in result
    assert result.count("blastqueryab($nodes1[$k],$nodes2[$l])") == 1
    assert result.count("blastqueryab($nodes2[$l],$nodes1[$k])") == 1
    assert "%blastquery=();" in result


def test_rejects_second_transformation():
    with pytest.raises(ValueError, match="already pair-parallel"):
        parallelize_source(parallelize_source(SERIAL_SOURCE))
