package Chemistry::Ring::Find;

$VERSION = 0.11;
# $Id: Find.pm,v 1.2 2004/06/18 00:38:00 ivan Exp $

=head1 NAME

Chemistry::Ring::Find

=head1 SYNOPSIS

    use Chemistry::Ring::Find ':all';

    # find the smallest ring containing $atom
    my $ring = find_ring($atom);

    # find all the rings containing $bond
    my @rings = find_ring($bond, all => 1);

    # see below for more options

=head1 DESCRIPTION

This module implements a simple breadth-first ring finding algorithm. It
does not find all the rings in the structure; it only finds the rings that 
include a starting atom or bond. Future versions may find all the rings in 
the molecule, or perhaps the Smallest Set of Smallest Rings.

=head1 FUNCTIONS

These functions may be exported explicitly, or all by using the :all tag, but
nothing is exported by default.

=over 

=cut


use strict;
use warnings;
use Chemistry::Ring;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(find_ring);
our %EXPORT_TAGS = ( all => \@EXPORT_OK ); 

our $DEBUG = 0;

=item find_ring($origin, %opts)

Find the smallest ring containg $origin, which may be either an atom or a bond.
Returns a Chemistry::Ring object. Options:

=over

=item all

If true, find all the rings containing $origin. If false, return the first ring
found. Defaults to false.

=item min

Only find rings with a the given minimum size. Defaults to zero.

=item max

Only find rings up to the given maximium size. Defaults to unlimited size.

=item size

Only find rings with this size. Same as setting min and max to the same size.
Default: unspecified.

=item exclude

An array reference containing a list of atoms that must NOT be present in the
ring. Defaults to the empty list.

=item mirror

If true, find each ring twice (forwards and backwards). Defaults to false.

=back

=cut

# $origin is an atom
# options: min, max, size, all, mirror, exclude
sub find_ring {
    my ($origin, %opts) = @_;
    my $min_size = $opts{min} || $opts{size} || 0;
    my $max_size = $opts{max} || $opts{size};
    my %paths;
    my %bond_paths;
    my @q;
    my @rings;
    my %used_end_nodes;
    my $required_bond;
    my %exclude;
    @exclude{ @{$opts{exclude} || []} } = ();

    if ($origin->isa("Chemistry::Bond")) {
        $required_bond = $origin;
        ($origin) = $origin->atoms;
    }
    @q = ($origin);
    $paths{$origin} = [$origin];
    $bond_paths{$origin} = [];
    # $path{$atom} means how to get to $atom from $origin

    my $a;
    while ($a = shift @q) {
        my $from = $paths{$a}[-2];
        print "at $a from $from\n" if $DEBUG;
        for my $bn ($a->bonds_neighbors($from)) {
            my $nei  = $bn->{to};
            my $bond = $bn->{bond};
            next if exists $exclude{$nei};
            print "  -> $nei\n" if $DEBUG;
            if ($paths{$nei}) {
                print "found a path collision... " if $DEBUG;
                # check to make sure that the ring really started at $origin
                # and the size is what was requested
                my $size = @{$paths{$nei}} + @{$paths{$a}} - 1;
                if($paths{$nei}[1] != $paths{$a}[1]
                    and $size >= $min_size
                    and !$max_size || $size <= $max_size)
                {
                    print "VALID\n" if $DEBUG;
                    my @atoms = (@{$paths{$a}}, reverse @{$paths{$nei}});
                    print "RING = ", print_path(\@atoms) if $DEBUG;
                    pop @atoms;
                    if ($used_end_nodes{$atoms[1]} and !$opts{mirror}) {
                        print "skipping redundant ring\n" if $DEBUG;
                        next; # don't want to find rings twice
                    }
                    my @bonds = (@{$bond_paths{$a}}, $bond,
                        reverse @{$bond_paths{$nei}});
                    if ($required_bond 
                        and not grep {$_ eq $required_bond} @bonds) {
                        print "does not include required bond\n" if $DEBUG;
                        next;
                    }
                    if (contains_ring(\@atoms, \@rings)) {
                        print "contains another ring\n" if $DEBUG;
                        next;
                    }
                    my $r = Chemistry::Ring->new;
                    $r->add_atom_np(@atoms);
                    $r->add_bond_np(@bonds);
                    return $r unless $opts{all};  # FOUND VALID RING
                    push @rings, $r;
                    $used_end_nodes{$atoms[-1]} = 1;
                    #@used_nodes{@atoms} = ();
                } else {
                    print "NOT VALID", 
                        print_path( [@{$paths{$a}}, 
                            reverse @{$paths{$nei->id}}]) if $DEBUG;
                }
            } else {
                if (!$max_size || @{$paths{$a}} < ($max_size / 2) + 0.1) {
                    push @q, $nei;
                    print "    pushing path\n" if $DEBUG;        
                    $paths{$nei} = [@{$paths{$a}}, $nei];
                    $bond_paths{$nei} = [@{$bond_paths{$a}}, $bond];
                    print print_path($paths{$nei}) if $DEBUG;
                } else {
                    print "path too long; " if $DEBUG;
                    print print_path($paths{$a}) if $DEBUG;
                    #path too long
                }
            }
        }
    }
    @rings;
}


sub print_path {
    my $p = shift;
    my $ret = "    PATH: ";
    for my $a (@$p) {
        $ret .=  "$a - ";
    }
    $ret .= "\n";
}

# contains_ring($atoms, $rings)
# returns true if one of the rings in the array ref $rings is a proper subset
# of the atom list in the array ref $atom.
sub contains_ring {
    my ($atoms, $rings) = @_;
    my %seen;
    @seen{@$atoms} = ();
    for my $ring (@$rings) {
        my $unique_atoms = $ring->atoms; 
        next if $unique_atoms >= @$atoms; # ring is same size or bigger
        # make sure that $ring has at least one atom not in $atoms
        for my $atom ($ring->atoms) {
            if (exists $seen{$atom}) {
                $unique_atoms--;
            } else {
                last; # there's at least one unique atom!
            }
        }
        return 1 unless $unique_atoms;
    }
    0;
}

1;

=back

=head1 VERSION

0.11

=head1 SEE ALSO

L<Chemistry::Ring>, L<Chemistry::Mol>, L<Chemistry::Atom>, L<Chemistry::Bond>
L<Math::VectorReal>.

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

