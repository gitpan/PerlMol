use Test::More;
use Chemistry::Mol;
use strict;

my @files;

BEGIN { 
    eval "use Chemistry::File::SMILES 0.40";
    @files = glob "t/*.pat" unless $@;
    plan tests => 1 + @files;
    use_ok('Chemistry::Pattern');
};


for my $file (@files) {
    open F, $file or die "couldn't open $file\n";   
    my ($patt_str, $options, $mol_str, @expected_matches) = map { /: ([^\n\r]*)/g } <F>;
    
    my ($mol, $patt);
    Chemistry::Atom->reset_id;
    $patt = Chemistry::Pattern->parse($patt_str, format => 'smiles');
    $mol = Chemistry::Mol->parse($mol_str, format => 'smiles');
    $patt->options(split " ", $options);

    my @matches;
    while ($patt->match($mol) ) {
        my @ret = $patt->atom_map;
        push @matches, "(@ret)";
    }
    push @matches, "()";

    is_deeply(\@matches, \@expected_matches, "$file: $mol_str =~ /$patt_str/");
}
