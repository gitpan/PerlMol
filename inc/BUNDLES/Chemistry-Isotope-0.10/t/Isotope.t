use Test::More;
use Data::Dumper;

#plan 'no_plan';
plan tests => 3;

my $TOL = 0.00001;

use Chemistry::Isotope qw(isotope_mass);

my $m = isotope_mass(235, 92);
is_float($m, 235.043923094753, "235U");

$m = isotope_mass(235, 6);
ok(! defined $m, "235C");

$m = isotope_mass(400, 180);
ok(! defined $m, "400, 180");

sub is_float {
    my ($got, $expected, $name) = @_;
    ok( abs($got-$expected) < $TOL, $name);
}
