use Test::More tests => 1;

use strict;
use warnings;
use Chemistry::File::XYZ;

my $mol = Chemistry::Mol->read("cyclobutene.xyz");
my $got = $mol->print(format => "xyz");
open F, "cyclobutene2.xyz" or die "Couldn't open cyclobutene2.xyz: $!\n";
my $expected = do { local $/; <F> };

is($got, $expected, "cyclobutene.xyz");
