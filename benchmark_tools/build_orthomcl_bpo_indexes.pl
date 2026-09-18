use strict;
use warnings;
use orthomcl_module;

die "Require absolute BPO and fresh absolute output directory\n" unless @ARGV == 2;
my ($bpo, $output) = @ARGV;
die "Require absolute BPO and fresh absolute output directory\n"
    unless $bpo =~ m{\A/} && $output =~ m{\A/} && -f $bpo && -s $bpo
        && !-e $output && !-l $output;
mkdir $output or die "Cannot create index directory: $!\n";
constructIDX_for_bpofile($bpo, "$output/all_bpo.idx");
constructSE_for_bpofile($bpo, "$output/all_bpo.se");
die "Native index output missing\n"
    unless -s "$output/all_bpo.idx" && -s "$output/all_bpo.se";
