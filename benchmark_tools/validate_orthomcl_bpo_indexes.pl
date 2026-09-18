use strict;
use warnings;
use Storable qw(retrieve);
use JSON::PP;

die "Require BPO, offset index and query-range index\n" unless @ARGV == 3;
my ($bpo, $idx, $se) = @ARGV;
my $offsets = retrieve($idx);
my $ranges = retrieve($se);
die "Invalid native index types\n" unless ref($offsets) eq 'ARRAY' && ref($ranges) eq 'HASH';
open(my $stream, '<:raw', $bpo) or die "Cannot read BPO: $!\n";
my ($count, $queries, $current, $start) = (0, 0, undef, undef);
my $finish_query = sub {
    my ($last) = @_;
    return unless defined $current;
    die "Wrong query range for $current\n" unless exists $ranges->{$current}
        && !ref($ranges->{$current}) && defined $ranges->{$current}
        && $ranges->{$current} eq "$start;$last";
    delete $ranges->{$current};
    $queries++;
};
while (1) {
    my $position = tell($stream);
    die "Wrong byte offset at entry $count\n" unless defined $offsets->[$count]
        && !ref($offsets->[$count]) && $offsets->[$count] =~ /\A(?:0|[1-9][0-9]*)\z/
        && $offsets->[$count] == $position;
    my $line = <$stream>;
    last unless defined $line;
    $line =~ s/\r?\n\z//;
    my @fields = split /;/, $line, -1;
    $count++;
    die "Malformed/nonsequential BPO at line $count\n"
        unless @fields == 8 && $fields[0] eq "$count" && length $fields[1];
    if (!defined $current || $fields[1] ne $current) {
        $finish_query->($count - 1);
        $current = $fields[1];
        $start = $count;
        die "Missing or repeated query block $current\n" unless exists $ranges->{$current};
    }
}
close($stream) or die $!;
die "Empty BPO\n" unless $count;
$finish_query->($count);
die "Extra query range entries\n" if keys %$ranges;
die "Wrong offset entry count\n" unless @$offsets == $count + 1;
print JSON::PP->new->canonical->encode({
    status => 'native_bpo_indexes_verified', records => 0 + $count,
    queries => 0 + $queries, offset_entries_including_eof => scalar @$offsets,
    bpo_bytes => 0 + (-s $bpo),
}), "\n";
