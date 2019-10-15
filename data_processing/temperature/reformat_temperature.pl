#! /usr/bin/perl

###############################################################################
# Writes temperature values into a format that is used to load into rhea.
# The temperature values correspond to vertices of the mesh from which
# coordinates were extracted.
#
# Usage: cat <input file> | perl reformat_temperature.pl > <input file>
###############################################################################

use strict;

my $line = 0;
while (<STDIN>) {
  if (!($_ =~ m/\s*(\d+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)/)) {
    die ("Line $line not matching\n");
  }
  printf ("%.6f\n", $5);
  ++$line;
}
