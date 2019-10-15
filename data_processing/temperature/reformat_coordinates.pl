#! /usr/bin/perl

###############################################################################
# Writes output of coordinates from rhea into a format that is provided as
# input to generate temperature values.
#
# Usage: cat <input file> | perl reformat_coordinates.pl > <output file>
###############################################################################

use strict;

my $line = 0;
while (<STDIN>) {
  if (!($_ =~ m/\s*([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s*/)) {
    die ("Line $line not matching\n");
  }
  printf ("%12d %+.3e %+.15e %+.15e\n", $line, $1, $2, $3);
  ++$line;
}
