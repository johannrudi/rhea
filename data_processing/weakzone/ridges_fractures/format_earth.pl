#!/bin/perl

# Usage: cat INFILE | perl format_earth.pl > OUTFILE

use strict;

my $line = 0;
LINE: while (<STDIN>) {
  next LINE if /^>/; # discard comments
  if (!($_ =~ m/\s*([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s*/)) {
    die ("ERROR: Line $line not matching\n");
  }
  printf("%+.6f %+.6f %+.6f\n", $1, $2, $3);
  ++$line;
}
