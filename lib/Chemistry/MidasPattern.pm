package Chemistry::MidasPattern;

$VERSION = '0.10';

# $Id: MidasPattern.pm,v 1.2 2004/07/03 01:05:22 itubert Exp $

=head1 NAME

Chemistry::MidasPattern - Select atoms in macromolecules

=head1 SYNOPSIS

    use Chemistry::MidasPattern;
    use Chemistry::File::PDB;

    # read a molecule
    my $mol = Chemistry::MacroMol->read("test.pdb");

    # define a pattern matching carbons alpha and beta
    # in all valine residues
    my $str  = ':VAL@CA,CB';
    my $patt = Chemistry::MidasPattern->new($str);

    # apply the pattern to the molecule
    $patt->match($mol);

    # extract the results
    for my $atom ($patt->atom_map) {
        printf "%s\t%s\n",  $atom->attr("pdb/residue_name"), $atom->name;
    }
    printf "FOUND %d atoms\n", scalar($patt->atom_map);

=head1 DESCRIPTION

This module partially implements a pattern matching engine for selecting atoms
in macromolecules by using Midas/Chimera patterns. See
L<http://www.cmpharm.ucsf.edu/~troyer/troff2html/midas/Midas-uh-3.html#sh-2.1>
for a detailed description of this language.

This module shares the same interface as L<Chemistry::Pattern>; to perform a
pattern matching operation on a molecule, follow these steps.

1) Create a pattern object, by parsing a string. Let's assume that the pattern
object is stored in $patt and that the molecule is $mol.

2) Execute the pattern on the molecule by calling $patt->match($mol).

3) If $patt->match() returns true, extract the "map" that relates the pattern
to the molecule by calling $patt->atom_map. These method returns a list of the
atoms in the molecule that are matched by the pattern. Thus $patt->atom_map(1)
would be analogous to the $1 special variable used for regular expresion
matching. The difference between Chemistry::Pattern and Perl regular
expressions is that atoms are always captured, and that each atom always uses
one "slot".

=head1 MIDAS ATOM SPECIFICATION LANGUAGE QUICK SUMMARY

The current implementation does not have the concept of a model, only of
residues and atoms.

What follows is not exactly a formal grammar specification, but it should give
a general idea: 

SELECTOR = ((:RESIDUE)*(@ATOM)*)*

The star here means "zero or more", and the parentheses are used to delimit the
effect of the star. The : and @ are used verbatim.

RESIDUE can be a name (e.g., LYS), a sequence number (e.g., 108), a range
(e.g., 1-10), or a comma-separated list of RESIDUEs (e.g. 1-10,6,LYS).

ATOM is an atom name, a serial number (this is a non-standard extension) or a
comma-separated list of ATOMs.

Names can have wildcards: * matches the whole name; ? matches one character;
and = matches zero or more characters. An @ATOM specification is asociated with
the closest preceding residue specification.

DISTANCE_SELECTOR = SELECTOR za< DISTANCE

Atoms within a certain distance of those that are matched by a selector can be
selected by using the za< operator, where DISTANCE is a number in Angstroms.

EXPR = ( SELECTOR | DISTANCE_SELECTOR ) (& (SELECTOR | DISTANCE_SELECTOR))*

The result of two or more selectors can be intersected using the & operator.

=head1 EXAMPLES

    :ARG                All arginine atoms
    :ARG@*              All arginine atoms
    @CA                 All alpha carbons
    :*@CA               All alpha carbons
    :ARG@CA             Arginine alpha carbons
    :VAL@C=             Valine carbons
    :VAL@C?             Valine carbons with two-letter names
    :ARG,VAL@CA         Arginine and valine alpha carbons
    :ARG:VAL@CA         All arginine atoms and valine alpha carbons
    :ARG@CA,CB          Arginine alpha and beta carbons
    :ARG@CA@CB          Arginine alpha and beta carbons
    :1-10               Atoms in residues 1 to 10
    :48-*               Atoms in residues 11 to the last one
    :30-40@CA & :ARG    Alpha carbons in residues 1-10 which are
                        also arginines.
    @123                Atom 123
    @123 za<5.0         Atoms within 5.0 Angstroms of atom 123
    @123 za>30.0        Atoms not within 30.0 Angstroms of atom 123
    @CA & @123 za<5.0   Alpha carbons within 5.0 Angstroms of atom 123

=cut


use strict;
use warnings;
use Chemistry::MacroMol;
use Carp;
use base "Chemistry::Pattern";

our $DEBUG = 0;

sub new {
    my ($class, $str) = @_;
    my $self = bless { 
        atom_map => [],
        options  => {},
        matched  => 0,
    }, ref $class || $class;
    $self->parse_midas_pattern($str);
    $self;
}

sub atom_map {
    my ($self) = @_;
    @{$self->{atom_map}};    
}

sub bond_map { () }

sub options {
    my $self = shift;
    $self->{options} = {@_};    
    $self;
}

# residue specifier: name | range | number | *
# atom specifier:    name | number | *
# names may have ? and = as wildcards (equivalent to . and .*)
sub parse_midas_pattern {
    my ($self, $string) = @_;
    my $tree = [];
    $self->{tree} = $tree;

    for my $str (split / +& +/, $string) {
        my $node = { type => 'normal' }; 
        if ($str =~ s/ +z(.?)([<>])\s*(\d+(\.\d+)?)//) {
            $node->{type} = 'zone_' . ($1 || 'r');
            $node->{zone_op} = $2;
            $node->{zone_val} = $3;
        }
        push @$tree, $node;
        my $conditions = [];
        $node->{conditions} = $conditions;

        $str = ":*$str" unless $str =~ /^:/; # add implicit "every residue"
        my (undef, @residues) = split ':', $str;# =~ /:([^:]+)/g;
        for my $residue (@residues) {
            my ($res_base, $atom_base) = split '@', $residue, 2;
            #print "res_base = $res_base\n";
            my @res = map { parse_res($_) } split ',', $res_base;
            my $res_node = { residue => \@res };
            push @$conditions, $res_node; 

            my (@atoms) = map { parse_atom($_) } split /[\@,]/, $atom_base||'';
            @atoms or @atoms = [all => '*'] ; # implicit "every atom"
            $res_node->{atoms} = \@atoms;
        }
    }
}

sub parse_res {
    my ($s) = @_;
    if ($s =~ /^\d+$/) {
        return [number => $s];
    } elsif ($s =~ /^[a-zA-Z]+$/) {
        return [name => $s];
    } elsif ($s =~ /^[a-zA-Z=?]+$/) {
        $s =~ s/\?/./g;
        $s =~ s/=/.*/g;
        return [regex => qr/$s/];
    } elsif ($s eq '*') {
        return [all => $s];
    } elsif ($s =~ /^(\d+|\*)-(\d+|\*)/) {
        return [range => $1, $2];
    } else {
        croak "invalid residue specification $s\n";
    }
}

sub parse_atom {
    my ($s) = @_;
    if ($s =~ /^\d+$/) {
        return [number => $s];
    } elsif ($s =~ /^[a-zA-Z0-9=?]+$/) {
        $s =~ s/\?/./g;
        $s =~ s/=/.*/g;
        return [regex => qr/^$s$/];
    } elsif ($s eq '*') {
        return [all => $s];
    } else {
        croak "invalid atom specification $s\n";
    }
}

sub match {
    my ($self, $mol) = @_;
    my $tree = $self->{tree};
    my @ret;
    my $first = 1;

    if ($self->{matched} and $self->{matched} eq $mol) {
        $self->{matched} = 0;
        return 0;
    }

    for my $term (@$tree) {
        print "term\n" if $DEBUG;        
        my @atoms;
        for my $residue_spec (@{$term->{conditions}}) {
            print "\tresidue_spec\n" if $DEBUG;        
            my @residues;
            
            # get the matching residues
            for my $cond (@{$residue_spec->{residue}}) {
                print "\t\tresidue condition\n" if $DEBUG;        
                my ($type, @params) = @$cond;
                if ($type eq 'all') {
                    push @residues, $mol->domains;
                } elsif ($type eq 'number') {
                    push @residues, $mol->domains(@params);
                } elsif ($type eq 'name') {
                    push @residues, grep {$_->type eq $params[0]} $mol->domains;
                } elsif ($type eq 'regex') {
                    push @residues, grep {$_->type =~ $params[0]} $mol->domains;
                } elsif ($type eq 'range') {
                    $params[0] = 1 if $params[0] eq '*';
                    $params[1] = $mol->domains if $params[1] eq '*';
                    push @residues, $mol->domains($params[0] .. $params[1]);
                }
                if ($self->{options}{unique}) {
                    my %seen;
                    @seen{@residues} = @residues;
                    @residues = values %seen;
                }
            }

            # get the matching atoms
            for my $res (@residues) {
                for my $cond (@{$residue_spec->{atoms}}) {
                    print "atom_cond\n" if $DEBUG;
                    my ($type, $param) = @$cond;
                    push @atoms, grep { 
                        my $match = 0;
                        if ($type eq 'regex') {
                            $match = $_->name =~ $param;
                        } elsif ($type eq 'all') {
                            $match = 1;
                        } elsif ($type eq 'number') {
                            $match = $_->attr("pdb/serial_number") == $param;
                        }
                    } $res->atoms;
                }
            }
        }
        my %seen;
        if ($self->{options}{unique}) {
            @seen{@atoms} = @atoms;
            @atoms = values %seen;
        }
        if ($term->{type} eq 'normal') {
            @seen{@atoms} = @atoms;
        } elsif ($term->{type} eq 'zone_r') {
            croak "zr operator not implemented\n";
        } elsif ($term->{type} eq 'zone_a') {
            my $d  = $term->{zone_val};
            my $gt = $term->{zone_op} eq '>';
            #print "zone_a ret(@ret)\n" unless $first;
            @atoms = grep {
                my $found = 0;
                for my $a (@atoms) {
                    $found = 1, last if ($a->coords - $_->coords)->length < $d;
                }
                $found xor $gt;
            } $first ? $mol->atoms : @ret;
            #print "zone_a atoms(@atoms)\n" unless $first;
            @seen{@atoms} = @atoms;
        }
        if ($first) {
            # keep all the atoms from this term
            @ret = @atoms;
            $first = 0;
        } else {
            # keep the intersection of the previous terms with the 
            # atoms from the current term
            my @tmp = @ret;
            @ret = ();
            for (@tmp) {
                push @ret, $_ if $seen{$_};
            }
        }
        last unless @ret; # don't bother evaluating the other terms
    }
    # optionally sort the result
    if ($self->{options}{sort}) {
        @ret =  map { $_->[0] } sort { $a->[1] <=> $b->[1] } 
            map { [$_, $_->attr("pdb/serial_number")] } @ret;
    }
    $self->{atom_map} = \@ret;
    $self->{matched} = $mol;
}

1;



=head1 CAVEATS

If a feature does not appear in any of the examples, it is probably not
implemented. For example, the zr< zone specifier, atom properties, Chimera
extensions such as chains, etc.

The zone specifiers (selection by distance) currently use a brute-force N^2
algorithm. You can optimize an & expression by putting the most unlikely
selectors first; for example

    :1-20 zr<10.0 & :38         atoms in residue 38 within 10 A of atoms
                                in residues 1-20 (slow)

    :38 & :1-20 zr<10.0         atoms in residue 38 within 10 A of atoms
                                in residues 1-20 (not so slow)
    
In the first case, the N^2 search measures the distance between every atom in
the molecule and every atom in residues 1-20, and then intersects the results
with the atom list of residue 28; the second case only measures the distance
between every atom in residue 38 with every atom in residues 1-20. The second
way is much, much faster for large systems.

Some day, a future version may implement a smarter algorithm...

=head1 VERSION

0.10

=head1 SEE ALSO

L<Chemistry::File::MidasPattern>, L<Chemistry::Pattern>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut



