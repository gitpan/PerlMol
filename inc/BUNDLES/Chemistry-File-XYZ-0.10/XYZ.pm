package Chemistry::File::XYZ;

$VERSION = '0.10';
# $Id: XYZ.pm,v 1.1.1.1 2004/06/30 22:03:30 itubert Exp $

use base qw(Chemistry::File);
use Chemistry::Mol;
use Carp;
use strict;
use warnings;

=head1 NAME

Chemistry::File::XYZ

=head1 SYNOPSIS

    use Chemistry::File::XYZ;

    # read a XYZ file
    my $mol = Chemistry::Mol->read("myfile.xyz");

    # write a XYZ file
    $mol->write("out.xyz");

=cut

=head1 DESCRIPTION

This module reads XYZ files. It automatically registers the 'xyz' format with
Chemistry::Mol, so that XYZ files may be identified and read by
Chemistry::Mol->read().

The XYZ format is not strictly defined and there are various versions floating
around; this module accepts the following:

First line: atom count (optional)

Second line: molecule name or comment (optional)

All other lines: symbol, x, y, and z coordinates separated by spaces, tabs, or
commas.

If the first line doesn't look like a number, the atom count is deduced from
the number of lines in the file. If the second line looks like it defines an
atom, it is assumed that there was no name or comment.

=cut

Chemistry::Mol->register_format(xyz => __PACKAGE__);

sub parse_string {
    my ($class, $s, %opts) = @_;

    my $mol_class  = $opts{mol_class}  || 'Chemistry::Mol';
    my $atom_class = $opts{atom_class} || $mol_class->atom_class;

    my @lines = split /(?:\n|\r\n?)/, $s;
    my $n_atoms;
    if ($lines[0] =~ /^\s*\d+\s*$/) {
        $n_atoms = shift @lines;
    }
    my $name = '';
    unless ($lines[0] =~ /^\s*[A-Z][a-z]?([,\s]+[eE\d.+-]+){3}/) {
        $name = shift @lines;
    }
    $n_atoms = @lines unless $n_atoms;

    my $mol = $mol_class->new(name => $name);

    my $i = 0;
    for my $line (@lines) {
        $i++;
        last if $i > $n_atoms;
        $line =~ s/^\s+//;
        my ($symbol, $x, $y, $z) = split /[\s,]+/, $line;
        $mol->new_atom(symbol => $symbol, coords => [$x, $y, $z]);
    }
    return $mol;
}

sub name_is {
    my ($class, $fname) = @_;
    $fname =~ /\.xyz$/i;
}

sub file_is {
    my ($class, $fname) = @_;
    $fname =~ /\.xyz$/i;
}

sub write_string {
    my ($class, $mol, %opts) = @_;

    %opts = (count => 1, name => 1, %opts);

    # header
    my $ret = '';
    $ret .= $mol->atoms . "\n" if $opts{count};
    if ($opts{name}) {
        my $name = defined $mol->name ? $mol->name : '';
        $name =~ s/\n.*//s;
        $ret .= "$name\n"
    }

    # body
    for my $atom ($mol->atoms) {
        $ret .= sprintf "%-2s %8.4f %8.4f %8.4f\n", $atom->symbol, 
            $atom->coords->array;
    }
    $ret;
}

1;

=head1 VERSION

0.10

=head1 SEE ALSO

L<Chemistry::Mol>, L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=cut

