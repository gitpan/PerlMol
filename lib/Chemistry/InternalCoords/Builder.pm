package Chemistry::InternalCoords::Builder;

$VERSION = '0.17';

use strict;
use warnings;
#use diagnostics;

use Chemistry::Bond::Find 'find_bonds';
use Chemistry::Canonicalize 'canonicalize';
use Chemistry::InternalCoords;
use List::Util qw(first reduce);
use base 'Exporter';

our @EXPORT_OK = qw(build_zmat);
our %EXPORT_TAGS = ( all => \@EXPORT_OK );


=head1 NAME

Chemistry::InternalCoords::Builder -  Build a Z-matrix from cartesian
coordinates

=head1 SYNOPSIS

    use Chemistry::InternalCoords::Builder 'build_zmat'; 

    # $mol is a Chemistry::Mol object
    build_zmat($mol);

=head1 DESCRIPTION

This module builds a Z-matrix from the cartesian coordinates of a molecule,
making sure that atoms are defined in a way that allows for efficient structure
optimizations and Monte Carlo sampling.

The algorithm tries to start at the center of the molecule and builds outward
in a breadth-first fashion. Improper dihedrals are used to ensure clean
rotation of groups without distortion. All distance and angle references use
real bonds and bond angles where possible (the exception being disconnected 
structures).

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=head1 FUNCTIONS

These functions may be exported, although nothing is exported by default.
To export all functions, use the ":all" tag.

=over

=item build_zmat($mol)

Build a Z-matrix from the cartesian coordinates of the molecule. Side effect
warning: This function modifies the molecule heavily! First, it finds the bonds
if there are no bonds defined already (for example, if the structure came from
and XYZ file with no bond information). Second, it canonicalizes the molecule,
as a means of finding the "topological center". Third, it builds the Z-matrix
using a breadth-first search. Fourth, it sorts the atoms in the molecule in the
order that they were defined in the Z-matrix.

=cut

sub build_zmat {
    my ($mol, %opts) = @_;

    find_bonds($mol) unless $mol->bonds;
    canonicalize($mol);

    my @atoms = sort { $b->attr("canon/class") <=> $a->attr("canon/class") }
        $mol->atoms;

    my $ats = [];
    for my $atom (@atoms) {
        next if $atom->attr("zmat/index");
        zmat_bfs($mol, $atom, $ats);
    }
    $mol->sort_atoms( sub {
        $_[0]->attr("zmat/index") <=> $_[1]->attr("zmat/index")
    });
}

sub zmat_bfs {
    my ($mol, $origin, $atoms) = @_;
    my @q = $origin;
    my $i = @$atoms;
    push @$atoms, $origin;
    $origin->attr("zmat/index", ++$i);
    add_ic($origin, $atoms);
    while (my $atom = shift @q) {
        #print "at $atom with $i\n";
        for my $nei (sorted_neighbors($atom)) {
            unless ($nei->attr("zmat/index")) {
                $nei->attr("zmat/index", ++$i);
                add_ic($nei, $atoms);
                push @q, $nei;
                push @$atoms, $nei;
            }
        }
    }
    $atoms;
}

# $atoms is the list of atoms that have been added so far
sub add_ic {
    my ($atom, $atoms) = @_;
    my $ic;
    my $n = $atom->attr("zmat/index");
    if ($n == 1) {
        $ic = Chemistry::InternalCoords->new($atom);
    } elsif ($n == 2) {
        $ic = Chemistry::InternalCoords->new($atom, find_length($atom, $atoms));
    } elsif ($n == 3) {
        $ic = Chemistry::InternalCoords->new(
            $atom, find_angle($atom, $atoms));
    } else {
        $ic = Chemistry::InternalCoords->new(
            $atom, find_dihedral($atom, $atoms));
    }
    $atom->internal_coords($ic);
}

# Choose a good length reference for $atom
sub find_length {
    my ($atom, $atoms) = @_;
    my $ref = first { $_->attr("zmat/index") } $atom->neighbors;
    unless ($ref) {
        $ref = ${ 
            reduce { $a->[0] < $b->[0] ? $a : $b } 
            map { [$atom->distance($_)] } 
            grep { $_ ne $atom } @$atoms;
        }[1];
    }
    ($ref, scalar $atom->distance($ref));
}

# Choose a good angle (and length) reference for $atom
sub find_angle {
    my ($atom, $atoms) = @_;
    my ($len_ref, $len_val) = find_length(@_);
    my $ang_ref = 
        first { 
            $_->attr("zmat/index") && $_ ne $len_ref && $_ ne $atom
        } $len_ref->neighbors($atom), @$atoms;
    my $ang_val = $atom->angle_deg($len_ref, $ang_ref);
    ($len_ref, $len_val, $ang_ref, $ang_val);
}

# Choose a good dihedral (and angle and length) reference for $atom
sub find_dihedral {
    my ($atom, $atoms) = @_;
    my ($len_ref, $len_val, $ang_ref, $ang_val) = find_angle(@_);
    my $dih_ref = 
        first { 
            $_->attr("zmat/index") && $_ ne $len_ref 
                && $_ ne $ang_ref && $_ ne $atom
        } $len_ref->neighbors($atom), $ang_ref->neighbors($len_ref), @$atoms;
    my $dih_val = $atom->dihedral_deg($len_ref, $ang_ref, $dih_ref);
    ($len_ref, $len_val, $ang_ref, $ang_val, $dih_ref, $dih_val);
}

sub sorted_neighbors {
    my ($atom) = @_;
    my @bn = $atom->neighbors;
    @bn = sort { 
        $b->attr("canon/class") <=> $a->attr("canon/class") 
    } @bn;
    @bn;
}

1;

=back

=head1 VERSION

0.17

=head1 CAVEATS

This version may not work properly for big molecules, because the
canonicalization step has a size limit.

=head1 TO DO

Some improvements for handling disconnected structures, such as making sure
that the intermolecular distance is short.

Allowing more control over how much the molecule will be modified: sort or not,
canonicalize or not...

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Atom>, L<Chemistry::InternalCoords>,
L<Chemistry::Bond::Find>, L<Chemistry::Canonicalize>, 
L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

