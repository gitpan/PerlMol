package Chemistry::File::Mopac;

$VERSION = '0.10';
# $Id: Mopac.pm,v 1.4 2004/06/01 20:21:46 itubert Exp $

use 5.006;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol 0.10;
use Chemistry::InternalCoords;
use List::Util 'first';
use Carp;

=head1 NAME

Chemistry::File::Mopac - MOPAC 6 input file reader

=head1 SYNOPSIS

    use Chemistry::File::Mopac;

    my $mol = Chemistry::Mol->read('file.mop');

=cut

=head1 DESCRIPTION

This module reads MOPAC 6 input files. It can handle both internal coordinates
and cartesian coordinates. It also extracts molecules from summary files,
defined as those files that match /SUMMARY OF/ in the third line. Perhaps a
future version will extract additional information such as the energy and
dipole from the summary file.

This module registers the C<mop> format with Chemistry::Mol. For detection
purposes, it assumes that filenames ending in .mop or .zt have the Mopac 
format, as well as files whose first line matches /am1|pm3|mndo|mdg|pdg/i 
(this may change in the future).

When the module reads an input file into $mol, it puts the keywords (usually
the first line of the file) in $mol->attr("mopac/keywords"), the comments
(usually everything else on the first three lines) in
$mol->attr("mopac/comments"), and the internal coordinates for each atom in
$atom->attr("internal_coords").  This part of the interface is not entirely
stable, as it is likely that the next version of Chemistry::Atom will include
explicit support for internal coordinates.

=cut

Chemistry::Mol->register_format("mop");

sub parse_string {
    my $class = shift;
    my $string = shift;
    my %opts = @_; 
    my $mol_class = $opts{mol_class} || "Chemistry::Mol";
    my $atom_class = $opts{atom_class} || "Chemistry::Atom";
    my $bond_class = $opts{bond_class} || "Chemistry::Bond";

    my $mol = $mol_class->new();
    my @lines = split "\n", $string;
    local $_;

    # do we have a summary file?
    if ($lines[2] =~ /SUMMARY OF/) { 
        # skip everything until FINAL GEOMETRY OBTAINED
        while (defined ($_ = shift @lines)) {
            last if /FINAL GEOMETRY OBTAINED/;
        }
    }

    my @header = splice @lines, 0, 3;
    my @keys;

    # crazy mopac header extension rules
    push @keys, $header[0];
    if ($keys[0] =~ /&/) {
        push @keys, $header[1];
        if ($keys[1] =~ /&/) {
            push @keys, $header[2];
        }
    } elsif ($keys[0] =~ /(\+\s|\s\+)/) {
        push @keys, $header[1];
        push @header, shift @lines;
        if ($keys[1] =~ /(\+\s|\s\+)/) {
            push @keys, $header[2];
            push @header, shift @lines;
        }
    }

    $mol->attr("mop/keywords" => join "\n", @keys);
    my $comment = join "\n", @header[@keys..$#header];
    $comment =~ s/\s+$//; 
    $mol->attr("mopac/comments" => $comment);
    $comment =~ s/\n/ /g;
    $mol->name($comment);
    
    # read coords
    my @coords;
    for (@lines) {
        # Sample line below
        # O    1.232010  1  128.812332  1  274.372818  1    3   2   1
        last if /^\s*$/; #blank line
        push @coords, [split];
    }
    
    # note: according to the MOPAC6 manual, triatomics must always use
    # internal coordinates; molecules that don't specify connectivity
    # (columns 7-9) use cartesian coordinates. I use column 8 for testing
    # because column 7 is sometimes used for the partial charge, adding to
    # the confusion. I assume that diatomics are always internal as well.
    if (@coords <= 3 or first {defined $_->[8]} @coords) { 
        read_internal_coords($mol, @coords);
    } else { # Cartesian coords
        read_cartesian_coords($mol, @coords);
    }

    return $mol;
}

sub read_internal_coords {
    my ($mol, @coords) = @_;

    my $i = 0; # atom index

    for my $coord (@coords) {
        $i++;
        my ($symbol, $len_val, $len_opt, $ang_val, $ang_opt,
            $dih_val, $dih_opt, $len_ref, $ang_ref, $dih_ref) = @$coord;

        # implicit links for the second and third atoms
        $len_ref ||= 1 if $i == 2;
        if ($i == 3) {
            $len_ref ||= 2;
            $ang_ref = $len_ref == 1 ? 2 : 1;
        }

        my $atom = $mol->new_atom(symbol => $symbol);
        my $ic = Chemistry::InternalCoords->new($atom,
            $len_ref, $len_val,
            $ang_ref, $ang_val, 
            $dih_ref, $dih_val);
        $atom->attr("internal_coords", $ic);
        $ic->add_cartesians;
    }
}

sub read_cartesian_coords {
    my ($mol, @coords) = @_;
    for my $coord (@coords) {
        my ($symbol, $x_val, $x_opt, $y_val, $y_opt,
            $z_val, $z_opt) = @$coord;
        my $atom = $mol->new_atom(symbol => $symbol,
            coords => [$x_val, $y_val, $z_val]);
    }
}


sub file_is {
    my ($class, $fname) = @_;
    
    return 1 if $fname =~ /\.(?:mop|zt)$/i;

    open F, $fname or croak "Could not open file $fname";
    
    my $line = <F>;
    close F;
    return 1 if $line =~ /am1|pm3|mndo|mdg|pdg/i;
    return 0;
}



1;

=head1 VERSION

0.10

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::File>, L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

