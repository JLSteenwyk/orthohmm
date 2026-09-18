use strict;
use warnings;
use Config;
use JSON::PP;
use POSIX;
use IO::Handle;
use Getopt::Long;
use File::Basename;
use Bio::SearchIO;
use Bio::SearchIO::blasttable;
use Bio::SeqIO;
use Bio::SeqIO::fasta;
use Storable;
use orthomcl_module;

my %mapped;
open(my $maps, '<', '/proc/self/maps') or die $!;
while (my $line = <$maps>) {
    my @fields = split /\s+/, $line, 6;
    next unless @fields == 6;
    chomp $fields[5];
    $mapped{$fields[5]} = 1 if $fields[5] =~ m{^/};
}
close($maps) or die $!;
print JSON::PP->new->canonical->pretty->encode({
    versions => {perl => "$^V", bioperl_searchio => "$Bio::SearchIO::VERSION",
                 storable => "$Storable::VERSION", orthomcl => "$orthomcl_module::VERSION"},
    config => {map {$_ => $Config{$_}} qw(archname version useithreads use64bitint use64bitall
                                      archlibexp privlibexp sitelibexp sitearchexp)},
    search_path => \@INC,
    loaded_modules => {map {$_ => $INC{$_}} keys %INC},
    mapped_files => [sort keys %mapped],
});
