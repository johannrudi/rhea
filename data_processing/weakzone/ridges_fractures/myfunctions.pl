# *********************************************
# converts radial coordinates (in degree) into Carthesian coordinates.
# Radial coordinates are -90 <= theta < 90 and -180 <= phi < 180
sub S2C {
  my ($x1, $x2, $x3);
  (my $phi, my $theta) = @_;
  # convert to standard spherical coords where 0 <= theta < 180
  $x1 = sin($theta/180*PI) * cos($phi/180*PI);
  $x2 = sin($theta/180*PI) * sin($phi/180*PI);
  $x3 = cos($theta/180*PI);
  return ($x1,$x2,$x3);
}

# converts Carthesian coordinates to radial theta, phi coordinates
sub C2S {
  my ($theta, $phi);
  my ($x1, $x2, $x3) = @_;
  $theta = atan2 (sqrt($x1* $x1 + $x2 * $x2), $x3);
  $phi = atan2($x2,$x1);
  return (180/PI*$phi,180/PI*$theta);
}


# derive norm of 3D vector
sub vec_norm {
  (my $x1, my $x2, my $x3) = @_;
  return sqrt( ($x1 * $x1 ) +
	       ($x2 * $x2 ) +
	       ($x3 * $x3 ));
}

# derive max iteratively;
sub max {
  my ($max, @vars) = @_;
  for (@vars) {
    $max = $_ if $_ > $max;
  }
  return $max;
}

# derive min iteratively;
sub min {
  my ($min, @vars) = @_;
  for (@vars) {
    $min = $_ if $_ < $min;
  }
  return $min;
}

# derive cross product between 2 vectors u and v
# input order: u1,u2,u3,v1,v2,v3.
sub cross_prod {
  my ($u1, $u2, $u3, $v1, $v2, $v3) = @_;
  my $n1 = $u2*$v3 - $u3*$v2;
  my $n2 = $u3*$v1 - $u1*$v3;
  my $n3 = $u1*$v2 - $u2*$v1;
  return ($n1, $n2, $n3);
}

# Derive the tangential vector of a spline curve
# at a given parameter value
sub curve_derivative {
  my ($spx1, $spx2, $spx3, $param) = @_;
  my ($lv, $rv);
  # 1st component
  $lv = $spx1->evaluate($param - EPS);
  $rv = $spx1->evaluate($param + EPS);
  my $dx1 = ($rv - $lv)/(2 * EPS);
  # 2nd component
  $lv = $spx2->evaluate($param - EPS);
  $rv = $spx2->evaluate($param + EPS);
  my $dx2 = ($rv - $lv)/(2 * EPS);
  # 3rd component
  $lv = $spx3->evaluate($param - EPS);
  $rv = $spx3->evaluate($param + EPS);
  my $dx3 = ($rv - $lv)/(2 * EPS);
  return ($dx1, $dx2, $dx3);
}

# weak zone function on reference domain
# Input:
# - depth: d
# - depth where dip angle changes: dipchange
# - shallow dip angle: dip1
# - intermediate dip angle: dip2
# Output: x and y values for weak zone
sub weakzone_fun {
  my ($d, $dipchange, $dip1, $dip2) = @_;
  $dip1 = $dip1/180*PI;
  $dip2 = $dip2/180*PI;
  my $sd = $dipchange * 0.25;   # smoothing param
  if ($d <= ($dipchange - $sd)) {
    return ( $d / tan($dip1) );
  } else {			# define smoothing polynomial
    my $a = (1/tan($dip2) - 1/tan($dip1)) / 4 / $sd;
    my $b = 1/tan($dip1) - 2 * $a * ($dipchange - $sd);
    my $c = ($dipchange - $sd) / tan($dip1) - $b * ($dipchange - $sd)
      - $a * ($dipchange - $sd) * ($dipchange - $sd);
    if ($d < $dipchange + $sd) {
      return ( $a * $d * $d + $b * $d + $c);
    } else {
      return ( $a * ($dipchange + $sd) * ($dipchange + $sd)
	       + $b * ($dipchange + $sd) + $c
	       + ($d - $dipchange - $sd) / tan($dip2));
    }
  }
}

# last expression returns true value
1;
