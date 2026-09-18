BEGIN {
    @INC = grep { !ref($_) && m{^/} } @INC;
}
use strict;
use warnings;

my $script = shift @ARGV;
die "Require an absolute native script path\n" unless defined $script && $script =~ m{^/} && -f $script;
$0 = $script;
local $@;
local $! = 0;
my $result = do $script;
die $@ if $@;
die "Cannot execute $script: $!\n" if !defined($result) && $!;
