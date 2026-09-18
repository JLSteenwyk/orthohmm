use strict;
use warnings;
use Bio::SeqIO;
use JSON::PP;
use Storable qw(retrieve);
use orthomcl_module;

my ($fasta, $blast, $output) = @ARGV;
die "Require FASTA, BLAST and fresh output directory\n" unless @ARGV == 3 && !-e $output;
mkdir $output or die "mkdir: $!\n";
my %lengths;
my $input = Bio::SeqIO->new(-file => $fasta, -format => 'fasta');
while (my $entry = $input->next_seq) {
    die "Duplicate or empty input\n" if exists $lengths{$entry->id} || !$entry->length;
    $lengths{$entry->id} = $entry->length;
}
die "Empty input\n" unless keys %lengths;
open(orthomcl_module::LOG, '>', "$output/native.log") or die $!;
blast_parse($blast, "$output/native.bpo", 1e-5, \%lengths,
            {format => 'blasttable', hsp => 1});
close(orthomcl_module::LOG) or die $!;
constructIDX_for_bpofile("$output/native.bpo", "$output/native.idx");
constructSE_for_bpofile("$output/native.bpo", "$output/native.se");
my $report = {
    offsets => retrieve("$output/native.idx"),
    query_ranges => retrieve("$output/native.se"),
    loaded_modules => {map {$_ => $INC{$_}} keys %INC},
    perl_version => "$^V", bioperl_searchio_version => "$Bio::SearchIO::VERSION",
    storable_version => "$Storable::VERSION", orthomcl_version => "$orthomcl_module::VERSION",
};
open(my $json, '>', "$output/native.json") or die $!;
print $json JSON::PP->new->canonical->pretty->encode($report);
close($json) or die $!;
