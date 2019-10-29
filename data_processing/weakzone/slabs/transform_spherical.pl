#!/bin/perl

# Usage: cat INFILE | perl transform_spherical.pl > OUTFILE

use strict;

# add `../perl_include` to the include path
use FindBin;
use lib "$FindBin::RealBin/../perl_include";

require "coordinates_functions.pl";

my ($line, $lon, $lat, $depth, $x, $y, $z);

$line = 0;
LINE: while (<STDIN>) {
  next LINE if /^>/; # discard comments
  if (!($_ =~ m/\s*([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s*/)) {
    die ("ERROR: Line $line not matching\n");
  }
  $lon = $1;
  $lat = $2;
  $depth = $3;
  ($x, $y, $z) = spherical_geo_to_cartesian($lon, $lat, $depth);
  printf("%+.6f %+.6f %+.6f\n", $x, $y, $z);
  ++$line;
}
