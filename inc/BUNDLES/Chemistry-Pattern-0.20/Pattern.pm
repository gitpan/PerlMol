package Chemistry::Pattern;
$VERSION = '0.20';
# $Id: Pattern.pm,v 1.16 2004/05/31 23:13:47 itubert Exp $

=head1 NAME

Chemistry::Pattern - Chemical substructure pattern matching

=head1 SYNOPSIS

    use Chemistry::Pattern;
    use Chemistry::Mol;
    use Chemistry::File::SMILES;

    # Create a pattern and a molecule from SMILES strings
    my $mol_str = "C1CCCC1C(Cl)=O";
    my $patt_str = "C(=O)Cl";
    my $mol = Chemistry::Mol->parse($mol_str, format => 'smiles');
    my $patt = Chemistry::Pattern->parse($patt_str, format => 'smiles');

    # try to match the pattern
    while ($patt->match($mol)) {
        @matched_atoms = $patt->atom_map;
        print "Matched: (@matched_atoms)\n";
        # should print something like "Matched: (a6 a8 a7)"
    }

=head1 DESCRIPTION

This module implements basic pattern matching for molecules.
The Chemistry::Pattern class is a subclass of Chemistry::Mol, so patterns
have all the properties of molecules and can come from reading the same
file formats. Of course there are certain formats (such as SMARTS)
that are exclusively used to describe patterns.

To perform a pattern matching operation on a molecule, follow these steps.

1) Create a pattern object, either by parsing a file or string, or by adding
atoms and bonds by hand by using Chemistry::Mol methods. Note that atoms and 
bonds in a pattern should be Chemistry::Pattern::Atom and 
Chemistry::Patern::Bond objects. Let's assume that the pattern object is 
stored in $patt and that the molecule is $mol.

2) Execute the pattern on the molecule by calling $patt->match($mol).

3) If $patt->match() returns true, extract the "map" that relates the pattern to
the molecule by calling $patt->atom_map or $patt->bond_map. These methods 
return a list of the atoms or bonds in the molecule that are matched by the 
corresponding atoms in the pattern. Thus $patt->atom_map(1) would be analogous
to the $1 special variable used for regular expresion matching. The difference
between Chemistry::Pattern and Perl regular expressions is that atoms and bonds
are always captured.

4) If more than one match for the molecule is desired, repeat from step (2) 
until match() returns false.

=cut

use 5.006;
use strict;
use Carp;
use base qw(Chemistry::Mol);
use Chemistry::Pattern::Atom;
use Chemistry::Pattern::Bond;

=head1 METHODS

=over 4

=item Chemistry::Pattern->new(name => value, ...)

Create a new empty pattern. This is just like the Chemistry::Mol constructor,
with one additional option: "options", which expects a hash reference (the
options themselves are described under the options() method).

=cut

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    $self->reset;
    $self->{options} = {overlap=>1, permute=>0};
    $self;
}

=item $pattern->options(option => value,...)

Available options:

=over

=item overlap

If true, matches may overlap. For example, the CC pattern could match twice
on propane if this option is true, but only once if it is false. This option
is true by default.

=item permute

Sometimes there is more than one way of matching the same set of pattern atoms
on the same set of molecule atoms. If true, return these "redundant" matches.
For example, the CC pattern could match ethane with two different permutations
(forwards and backwards). This option is false by default.

=back

=cut

sub options {
    my $self = shift;
    if (@_ == 1) {
        $self->{options} = {%{$self->{options}}, %{$_[0]}};
    } else {
        $self->{options} = {%{$self->{options}}, @_};
    }
}

=item $mol->new_atom(name => value, ...)

Shorthand for $mol->add_atom(Chemistry::Atom->new(name => value, ...));

=cut

sub new_atom {
    my $self = shift;
    $self->add_atom(Chemistry::Pattern::Atom->new(@_));
}


=item $mol->new_bond(name => value, ...)

Shorthand for $mol->add_bond(Chemistry::Bond->new(name => value, ...));

=cut

sub new_bond {
    my $self = shift;
    $self->add_bond(Chemistry::Pattern::Bond->new(@_));
}

our $DEBUG = 0;

sub reset {
    my ($self, @mols) = @_;
    print "Resetting to (@mols)\n" if $DEBUG;
    $self->{pending_mols} = [@mols];
    $self->{pending_atoms} = [];
    $self->{already_matched} = {};
    $self->{current_atom} = $self->next_atom;
    $self->{paint_tab} = {};
    $self->match_local_init($self->{current_atom});
    $self->{flat} = $self->flatten;
}

sub already_matched {
    my ($self, @objs) = @_;
    my @ids = map {$_->id} @objs;
    my $unsorted_key = join " ", @ids;
    my $key;
    if ($self->{options}{permute}) {
        $key = $unsorted_key;
    } else {
        $key = join " ", sort @ids;
    }

    if ($self->{already_matched}{$key}) {
        print "already matched $unsorted_key\n" if $DEBUG;
        return 1;
    } else {
        $self->{already_matched}{$key} = 1;
        print "first match of $key\n" if $DEBUG;
        return 0;
    }
}

sub next_atom {
    my $self = shift;
    my $atom;
    print "next_atom\n" if $DEBUG;
    if (@{$self->{pending_atoms}}) {
        $atom = shift @{$self->{pending_atoms}};
        print "\tatom $atom\n" if $DEBUG;
    } elsif (@{$self->{pending_mols}}) {
        my $mol = shift @{$self->{pending_mols}};
        print "\tmol $mol\n" if $DEBUG;
        $self->map_to($mol);
        $self->{pending_atoms} = [$mol->atoms];
        $atom = shift @{$self->{pending_atoms}};
        print "\tatom $atom\n" if $DEBUG;
    }
    $atom;
}

Chemistry::Obj::accessor "map_to";

=item $pattern->atom_map

Returns the list of atoms that matched the last time $pattern->match was called.

=cut

sub atom_map {
    my $self = shift;
    my @atoms = map { $_->map_to } $self->atoms(@_);
    @atoms == 1 ? $atoms[0] : @atoms;
}

=item $pattern->bond_map

Returns the list of bonds that matched the last time $pattern->match was called.

=cut

sub bond_map {
    my $self = shift;
    my @bonds = map { $_->map_to } $self->bonds(@_);
    @bonds == 1 ? $bonds[0] : @bonds;
}

=item $pattern->match($mol)

Returns true if the pattern matches the molecule. If called again for the 
same molecule, continues matching where it left off (in a way similar to global
regular expressions under scalar context). When there are no matches left,
returns false.

To find out which atoms and bonds matched, use the atom_map and bond_map
methods.

=cut

# testing a different implementation
sub match {
    my ($self, $mol) = @_;
    print "match $self $mol\n" if $DEBUG;
    if (defined($mol) and $self->map_to ne $mol) { #TODO
        $self->reset($mol);
    }
    my $match = $self->match_next;
    print "returning match: '$match'\n" if $DEBUG;
    $match;
}

sub match_next {
    my $self = shift;
    my $match;
    print "match_next\n" if $DEBUG;

    $self->match_local_init($self->{current_atom}) 
        unless $self->{options}{overlap};
    while (1) {
        $match = $self->match_local_next;
        if ($match) {
            if ($self->already_matched($self->atom_map, $self->bond_map)) {
                $match = 0;
                next; # already matched this; try again
            } else {
                unless ($self->{options}{overlap}) {
                    $self->next_atom;
                }
                last; # matched!
            }
        } else { 
            my $atom = $self->next_atom;
            last unless $atom;
            $self->match_local_init($atom);
        }
    }
    $match;
}


sub match_local_init {
    my ($patt, $atom) = @_;
    my $state = 1;
    my @stack = {
        from_where_bond_i => [0], 
        from_what_bond_i => [0],
        from_where => [$atom], 
        from_what => [$patt->atoms(1)],
    };
    my $mol = $patt->map_to;

    print "match_local_init(",$patt->atoms(1),", $atom)\n" if $DEBUG;
    for ($patt->atoms, $patt->bonds) { $_->map_to(undef) }
    if ($mol and $patt->{options}{overlap}) {
        $patt->{paint_tab} = {};
    }
    $patt->{stack} = \@stack;
    $patt->{state} = $state;

    # new stuff
    $patt->{stack2} = [[0]];
    $patt->{backtrack} = 0;
    $patt->{current_atom} = $atom;
}


sub match_local_next {
    my ($patt) = @_;
    my $stack = $patt->{stack2};
    my $paint_tab = $patt->{paint_tab};
    my $match = 0;
    my $backtrack = $patt->{backtrack};
    my $where = $patt->{current_atom};
    my $flat = $patt->{flat};
    my $mol = $patt->map_to;

    print "match_local_next($where)\n" if $DEBUG;

    while (1) {
        my ($idx) = @{$stack->[-1]};
        my $pos = @$stack - 1;
        my $what = $flat->[$pos];
        if($DEBUG) { print "[@$_] " for @$stack; print "\n" }

        if ($backtrack) {
            print "\tbacktracking...\n" if $DEBUG;
            if (@$stack <= 2) {
                print "\tnowhere to backtrack\n" if $DEBUG;
                last;
            }
            my $to;
            pop @$stack;
            (undef, $what, $to) = @{$stack->[-1]};
            $backtrack = 0;
            if ($what) {
                print "\tcleaning $what\n" if $DEBUG;
                $paint_tab->{$what->map_to} = 0;
                $what->map_to(undef);
            }
            if ($to) {
                print "\tcleaning $to\n" if $DEBUG;
                $paint_tab->{$to->map_to} = 0;
                $to->map_to(undef);
            }
            next;
        } elsif ($pos >= @$flat) { # nothing else to match
            print "\tmatched!\n" if $DEBUG;
            $backtrack = 1;
            $match = 1;
            last;
        } elsif (@$stack == 1) { # match first atom (anchored)
            if (!$paint_tab->{$where} and $what->test($where)) {
                print "\tmatched initial atom $what on $where\n" if $DEBUG;
                $paint_tab->{$where} = 1;
                $what->map_to($where);
                push @$stack, [0];
                next;
            } else {
                last;
            }
        } elsif ($what->isa("Chemistry::Pattern::Atom")) {
            # match unanchored atom
            if($idx >= $mol->atoms) {
                print "\tno more atoms at $mol\n" if $DEBUG;
                $backtrack = 1, next;
            }
            my $next = $mol->atoms($idx+1);
            $stack->[-1] = [++$idx];
            if ($paint_tab->{$next}) {
                print "\t$next already visited\n" if $DEBUG;
                next;
            }
            if ($what->test($next)) {
                print "\tmatched unanchored atom $what on $next\n" if $DEBUG;
                $paint_tab->{$next} = 1;
                $what->map_to($next);
                $stack->[-1] = [$idx, undef, $what];
                push @$stack, [0];
            } else {
                print "\tunanchored atom $what didn't match $next\n" if $DEBUG;
            }
            next;
        } elsif ($what->isa("Chemistry::Pattern::Bond")) { # match bond
            my ($a1, $a2) = $what->atoms;
            my ($from, $to) = $a1->map_to ? ($a1, $a2) : ($a2, $a1);
            $where = $a1->map_to || $a2->map_to;
            print "\tfrom $from to $to\n" if $DEBUG;
            my @bn = $where->bonds_neighbors;
            if ($idx >= @bn) {
                print "\tno more bonds at $where\n" if $DEBUG;
                $backtrack = 1, next;
            }
            my $nei  = $bn[$idx]{to};
            my $bond = $bn[$idx]{bond};
            $stack->[-1] = [++$idx];

            if ($paint_tab->{$bond}) {
                print "\t$bond already visited\n" if $DEBUG;
                next;
            }
            if ($a1->map_to and $a2->map_to) { # ring closure bond
                if ($what->test($bond)) {
                    print "\tmatched ring bond $what on $bond\n" if $DEBUG;
                    if ($to->map_to eq $nei) {
                        print "\tring closed $to on $nei\n" if $DEBUG;
                        $paint_tab->{$bond} = 1;
                        $what->map_to($bond);
                        $stack->[-1] = [$idx, $what];
                        push @$stack, [0];
                        next;
                    } else {
                        print "\tring didn't close $to on $nei\n" if $DEBUG;
                    }
                }
            } else {
                if ($paint_tab->{$nei}) {
                    print "\t$nei already visited\n" if $DEBUG;
                    next;
                }
                if ($what->test($bond)) {
                    print "\tmatched bond $what on $bond\n" if $DEBUG;
                    if ($to->test($nei)) {
                        print "\tmatched atom $to on $nei\n" if $DEBUG;
                        $paint_tab->{$nei} = 1;
                        $to->map_to($nei);
                        $paint_tab->{$bond} = 1;
                        $what->map_to($bond);
                        $stack->[-1] = [$idx, $what, $to];
                        push @$stack, [0];
                        next;
                    } else {
                        print "\tatom $to didn't match $nei\n" if $DEBUG;
                    }
                } else {
                    print "\tbond $what didn't match $bond\n" if $DEBUG;
                }
            }
            next;
        } else {
            die "shouldn't be here!";
        }
    }
    $patt->{backtrack} = $backtrack;
    $patt->{stack2} = $stack;
    $match;
}

sub flatten {
    my ($patt) = @_;
    my $visited = {};
    my $list = [];
    for my $atom ($patt->atoms) {
        next if $visited->{$atom};
        push @$list, $atom;
        _flatten($patt, $atom, $visited, $list); 
    }
    print "flattened list: (@$list)\n" if $DEBUG;
    $list;
}

sub _flatten {
    my ($patt, $atom, $visited, $list) = @_;

    $visited->{$atom}  = 1;
    for my $bn ($atom->bonds_neighbors) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond}  = 1;
        push @$list, $bond;
        unless ($visited->{$nei}) {
            _flatten($patt, $nei, $visited, $list);
        }
    }
}

1;

=back

=head1 VERSION

0.20

=head1 SEE ALSO

L<Chemistry::Pattern::Atom>, L<Chemistry::Pattern::Bond>, L<Chemistry::Mol>,
L<Chemistry::File>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

