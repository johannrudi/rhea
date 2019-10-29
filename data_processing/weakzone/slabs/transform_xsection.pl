#!/bin/perl

###############################################################################
# Usage: 
#   cat <input file> | \
#   perl transform_xsection.pl <resolution in km> > <output file>
###############################################################################

use strict;
use List::Util qw[min max];

# add `../perl_include` to the include path
use FindBin;
use lib "$FindBin::RealBin/../perl_include";

require "coordinates_functions.pl";

# process script arguments
if ($#ARGV != 0) {
    die "Usage: $0 <resolution in km>";
}
my $resolution_km = $ARGV[0];

# set parameters
my $pi = 3.14159265359;
my $earth_radius_km = 6371;
my $resolution = $resolution_km/$earth_radius_km;
my $lon_min = -$pi/8.0 * 82.0/16.0;
my $lon_offset = $pi/2.0 + $lon_min;
my $lat_min = -$pi/8.0 * 1.0/16.0;
my $lat_max = +$pi/8.0 * 1.0/16.0;
my $n_lat = 1 + max(0, int(($lat_max - $lat_min)/$resolution));

# declare variables
my ($line, $lon_deg, $depth_km, $lat, $x, $z, $x_prev, $z_prev,
    $dx, $dz, $hx, $hz, $dist, $n_interpolate, $xi, $yi, $zi);

$line = 0;
LINE: while (<STDIN>) {
  # read line
  next LINE if /^>/; # discard comments
  if (!($_ =~ m/\s*([-+\d\.e]+)\s+([-+\d\.e]+)\s*/)) {
    die ("ERROR: Line $line not matching\n");
  }
  $lon_deg = $1 + 180.0*$lon_offset/$pi;
  $depth_km = (1.0 - $2)*$earth_radius_km;

  # transform (lon,depth) coordinates to (x,z)
  ($x, $z) = polar_geo_to_cartesian($lon_deg, $depth_km);

  # set (x,z)-interpolation parameters
  $n_interpolate = 0;
  $hx = 0.0;
  $hz = 0.0;
  if (0 < $line) {
    $dx = $x - $x_prev;
    $dz = $z - $z_prev;
    $dist = sqrt($dx*$dx + $dz*$dz);
    if ($resolution < $dist) {
      $n_interpolate = max(0, int($dist/$resolution - 1));
    }
    $hx = $dx/(1.0 + $n_interpolate);
    $hz = $dz/(1.0 + $n_interpolate);
  }

 #printf("###DEV### line %i #interp %i #lat %i\n", $line, $n_interpolate, $n_lat);
  for (my $i = $n_interpolate; 0 <= $i; --$i) {
    # interpolate (x,z)-coordinates
    $xi = $x - $hx*$i;
    $zi = $z - $hz*$i;
    # interpolate y-coordinates
    for (my $j = 0; $j <= $n_lat; ++$j) {
      $lat = $lat_min + ($lat_max - $lat_min)*$j/$n_lat;
      $yi = sin($lat);
      printf("%+.12f %+.12f %+.12f\n", $xi, $yi, $zi);
    }
  }

  # proceed to next line
  $x_prev = $x;
  $z_prev = $z;
  ++$line;
}
