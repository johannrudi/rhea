###############################################################################
# Author:             Johann Rudi
###############################################################################

use Math::Trig;
use Carp::Assert;
use Machine::Epsilon;

# Converts spherical coordinates (lon_deg,lat_deg,depth_km) into Cartesian
# coordinates (x,y,z), where
# 
#   lon_deg  = positive [0,180] when east of prime meridian, 
#              negative [-180,0] when west
#   lat_deg  = positive [0,90] when north of equator, 
#              negative [-90,0] when south
#   depth_km = 0 at surface, 6371 at earth center
sub spherical_geo_to_cartesian {
  my ($rho, $theta, $phi, $x, $y, $z);

  # get & check input
  (my $lon_deg, my $lat_deg, my $depth_km) = @_;
  assert(-180 <= $lon_deg && $lon_deg <= 180) if DEBUG;
  assert(-90 <= $lat_deg && $lat_deg <= 90) if DEBUG;
  assert(0 <= $depth_km && $depth_km <= 6371) if DEBUG;

  # nondimensionalize input
  $rho = (6371.0 - $depth_km)/6371.0;
  $theta = pi * ($lon_deg + 180.0)/180.0;
  $phi = pi * (90.0 - $lat_deg)/180.0;
 #printf("###DEV### %+.8e %+.16e %+.8e\n", 0, pi*2, pi);
 #printf("###DEV### %+.8e %+.16e %+.8e\n", $rho, $theta, $phi);
  assert(0 <= $rho && $rho <= 1) if DEBUG;
  assert(0 <= $theta && $theta < pi*2 + 1000*machine_epsilon()) if DEBUG;
  assert(0 <= $phi && $phi < pi + 1000*machine_epsilon()) if DEBUG;

  # transform to Cartesian coordiantes
  #$x = $rho * cos($theta) * sin($phi);
  #$y = $rho * sin($theta) * sin($phi);
  #$z = $rho               * cos($phi);
  ($x, $y, $z) = Math::Trig::spherical_to_cartesian($rho, $theta, $phi);
  assert(-1 <= $x && $x <= 1) if DEBUG;
  assert(-1 <= $y && $y <= 1) if DEBUG;
  assert(-1 <= $z && $z <= 1) if DEBUG;

  # return Cartesian coordinates
  return ($x, $y, $z);
}

# Converts Cartesian coordinates (x,y,z) into spherical coordinates
# (lon_deg,lat_deg,depth_km).  
# Inverse function to `spherical_geo_to_cartesian()` above.
sub cartesian_to_spherical_geo {
  my ($rho, $theta, $phi, $lon_deg, $lat_deg, $depth_km);

  # get & check input
  (my $x, my $y, my $z) = @_;
  assert(-1 <= $x && $x <= 1) if DEBUG;
  assert(-1 <= $y && $y <= 1) if DEBUG;
  assert(-1 <= $z && $z <= 1) if DEBUG;

  # transform to spherical coordiantes
  #$rho = sqrt($x*$x + $y*$y + $z*$z);
  #$theta = atan2($y, $x);
  #$phi = acos($z / $rho);
  ($rho, $theta, $phi) = Math::Trig::cartesian_to_spherical($x, $y, $z);
  assert(0 <= $rho && $rho <= 1) if DEBUG;
  assert(0 <= $theta && $theta < pi*2) if DEBUG;
  assert(0 <= $phi && $phi < pi) if DEBUG;

  # dimensionalize according to geophysics convention
  $lon_deg = 180 * ($theta/pi - 1);
  $lat_deg = 180 * (0.5 - $phi/pi);
  $depth_km = 6371 * (1 - $rho);
  assert(-180 <= $lon_deg && $lon_deg <= 180) if DEBUG;
  assert(-90 <= $lat_deg && $lat_deg <= 90) if DEBUG;
  assert(0 <= $depth_km && $depth_km <= 6371) if DEBUG;

  # return spherical coordinates
  return ($lon_deg, $lat_deg, $depth_km);
}

# Converts polar coordinates (lon_deg,radius) into Cartesian coordinates (x,y).
sub polar_geo_to_cartesian {
  my ($rho, $theta, $x, $y, $z);

  # get & check input
  (my $lon_deg, my $depth_km) = @_;
 #assert(0 <= $lon_deg && $lon_deg <= 360) if DEBUG;
  assert(0-6371*0.001 <= $depth_km && $depth_km <= 6371*1.001) if DEBUG;

  # nondimensionalize input
  $rho = (6371.0 - $depth_km)/6371.0;
  $theta = 2.0*pi * $lon_deg/360.0;
  assert(0 <= $rho && $rho <= 1.001) if DEBUG;
 #assert(0 <= $theta && $theta < pi*2 + 1000*machine_epsilon()) if DEBUG;

  # transform to Cartesian coordiantes
  #$x = $rho * cos($theta);
  #$y = $rho * sin($theta);
  ($x, $y, $z) = Math::Trig::cylindrical_to_cartesian($rho, $theta, 0.0);
  assert(-1*1.003 <= $x && $x <= 1*1.003) if DEBUG;
  assert(-1*1.003 <= $y && $y <= 1*1.003) if DEBUG;

  # return Cartesian coordinates
  return (-$x, $y);
}

# last expression returns true value
1;
