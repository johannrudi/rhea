#! /usr/bin/perl

###############################################################################
# Lists depths of vertices (in km) of a spherical mesh from rhea, given a
# uniform level of mesh refinement.
#
# Usage: perl list_depths.pl <refinement level> <max depth in km>
###############################################################################

use strict;
use POSIX;

my $earth = 6371;
my $inrad = 0.55;
my $outrad = 1.0;

if ($#ARGV != 1) {
    die "Usage: $0 <refinement level> <max depth in km>";
}
my $level = $ARGV[0];
my $maxdepth = $ARGV[1];

my $tpowl = 2 ** $level;
my $n = $tpowl + 1;
my ($i, $rad, $z);

for ($i = $tpowl; $i >= 0; --$i) {
  $z = $earth * (1 -
      $inrad * $inrad / $outrad
      * pow ($outrad / $inrad, (1 + $i / $tpowl)));
  #$z = floor ($z + .5);
  if ($z <= $maxdepth) {
    printf "%04.3f\n", $z;
  }
  else {
    last;
  }
}
