#!/bin/perl

###############################################################################
# Usage: 
#   cat <input file> | \
#   perl generate_xsection_coordinates.pl <resolution in km> > <output file>
###############################################################################

use strict;
use List::Util qw[min max];

# add `../perl_include` to the include path
use FindBin;
use lib "$FindBin::RealBin/../perl_include";

require "coordinates_functions.pl";

# process script arguments
if ($#ARGV != 1) {
    die "Usage: $0 <resolution in km> <depth in km>";
}
my $resolution_km = $ARGV[0];
my $depth_max_km = $ARGV[1];

# set parameters
my $pi = 3.14159265359;
my $earth_radius_km = 6371;
my $resolution = $resolution_km/$earth_radius_km;
my $lon_min = -$pi/8.0 * 82.0/16.0;
my $lon_offset = $pi/2.0 + $lon_min;
my $lat_min = -$pi/8.0 * 1.0/16.0;
my $lat_max = +$pi/8.0 * 1.0/16.0;
my $n_lat = 1 + max(0, int(($lat_max - $lat_min)/$resolution));
my $n_depth = 1 + max(0, int(($depth_max_km/$earth_radius_km)/$resolution));

# declare variables
my ($line, $lon_deg, $depth_km, $lat, $xi, $yi, $zi);

$line = 0;
LINE: while (<STDIN>) {
  # read line
  next LINE if /^>/; # discard comments
  if (!($_ =~ m/\s*([-+\d\.e]+)\s*/)) {
    die ("ERROR: Line $line not matching\n");
  }
  $lon_deg = $1 + 180.0*$lon_offset/$pi;

 #printf("###DEV### line %i #depth %i #lat %i\n", $line, $n_depth, $n_lat);
  for (my $i = 0; $i <= $n_depth; ++$i) {
    # interpolate (x,z)-coordinates
    $depth_km = $depth_max_km*$i/$n_depth;
    ($xi, $zi) = polar_geo_to_cartesian($lon_deg, $depth_km);
    # interpolate y-coordinates
    for (my $j = 0; $j <= $n_lat; ++$j) {
      $lat = $lat_min + ($lat_max - $lat_min)*$j/$n_lat;
      $yi = sin($lat);
      printf("%+.12f %+.12f %+.12f\n", $xi, $yi, $zi);
    }
  }

  # proceed to next line
  ++$line;
}
