
# reads in data from Trench files
# Input: Datafile
#
# Output:
# nzones: Number of zones
# names: Names of zones
# lnum: line numbers of first point in each zone
# tvec, pvec: theta and phi values of surface line data
# off: Lateral surface shift for each point
# sdip: shallow dipping angle
# dep: depth in which dipping angle changes
# idip: intermediate dipping angle
# stop: overall depth of weak zone

sub readtrench {
  my ($inputfile) = @_;
  print "Inputfile: ".$inputfile."\n";
  open (IN, $inputfile) or die "Cannot open $inputfile\n";

  my (@off0, @off1, @sdip0, @sdip1, @dep0, @dep1, @idip0, @idip1,
      @stop0, @stop1, @reduction);
  my (@off, @sdip, @dep, @idip, @stop);
  my (@tvec, @pvec, @lnum, @names);

  my ($nzones, $npts);
  my ($pts, $offd, $sdipd, $depd, $idipd, $stopd);

  my $zone = 0;
  my $shift = 0;
  my ($dd1, $dd2, $ddd);
  $line = 0;

  while (<IN>) {
    if ($_ =~ m/\s?([-+\d\.eE]+)\s+([-+\d\.eE]+)\s?/) {
      $pvec[$line] = $1 - $shift;
      $tvec[$line] = $2;
      $line++;
    } elsif ($_=~ m/\s?.*sL\s?(\w+).*\s+OFF=(\d+),(\d+)\s?S_DIP=(\d+),(\d+)\s?DEPTH=(\d+),(\d+)\s?I_DIP=(\d+),(\d+)\s?V_REDUC=([-+\d\.eE]+)\s?/) {
      $lnum[$zone] = $line;
      $names[$zone] = $1;
      $off0[$zone] = $2;
      $off1[$zone] = $3;
      $sdip0[$zone] = $4;
      $sdip1[$zone] = $5;
      $dep0[$zone] = $6;
      $dep1[$zone] = $7;
      $idip0[$zone] = $8;
      $idip1[$zone] = $9;
      $stop0[$zone] = 120;
      $stop1[$zone] = 120;
      $reduction[$zone] = $10;
      print "Reading data for ".$1."\n";
      $zone++;
    } elsif ($_=~ m/\s?.*RI\s?(\w+).*\s+MAX_DEPTH=(\d+),(\d+)\s?/) {
      $lnum[$zone] = $line;
      $names[$zone] = $1;
      $off0[$zone] = 0;
      $off1[$zone] = 0;
      $sdip0[$zone] = 90;
      $sdip1[$zone] = 90;
      $dep0[$zone] = 50;
      $dep1[$zone] = 50;
      $idip0[$zone] = 90;
      $idip1[$zone] = 90;
      $stop0[$zone] = $2;
      $stop1[$zone] = $3;
      $reduction[$zone] = 0;
      print "Reading data for ridge ".$1."\n";
      $zone++;
    } elsif ($_=~ m/\s?.*FZ\s?(\w+).*\s+MAX_DEPTH=(\d+),(\d+)\s?/) {
      $lnum[$zone] = $line;
      $names[$zone] = $1;
      $off0[$zone] = 0;
      $off1[$zone] = 0;
      $sdip0[$zone] = 90;
      $sdip1[$zone] = 90;
      $dep0[$zone] = 50;
      $dep1[$zone] = 50;
      $idip0[$zone] = 90;
      $idip1[$zone] = 90;
      $stop0[$zone] = $2;
      $stop1[$zone] = $3;
      $reduction[$zone] = 0;
      print "Reading data for fracture ".$1."\n";
      $zone++;
    } else {
      die "Cannot read line $line";
    }
  }
  $lnum[$zone] = $line;
  $nzones = $zone;
  close($inputfile);
# interpolate offsets and dipping angles
  for ($zone = 0; $zone < $nzones; ++$zone) {
      $npts = $lnum[$zone + 1] - $lnum[$zone];
      print "Data for ".$names[$zone]." contains ".$npts.
	  " points.\n";
      my $offset = $lnum[$zone];
      for ($pts = 0; $pts < $npts; ++$pts) {
	  $off[$offset + $pts] = $off0[$zone] + ($pts/($npts-1)) *
	      ($off1[$zone] - $off0[$zone]);
	  $sdip[$offset + $pts] = $sdip0[$zone] + ($pts/($npts-1)) *
	      ($sdip1[$zone] - $sdip0[$zone]);
	  $dep[$offset + $pts] = $dep0[$zone] + ($pts/($npts-1)) *
	      ($dep1[$zone] - $dep0[$zone]);
	  $idip[$offset + $pts] = $idip0[$zone] + ($pts/($npts-1)) *
	      ($idip1[$zone] - $idip0[$zone]);
	  $stop[$offset + $pts] = $stop0[$zone] + ($pts/($npts-1)) *
	      ($stop1[$zone] - $stop0[$zone]);
      }
  }
# return references to arrays
  return ($nzones, \@names, \@lnum, \@pvec, \@tvec, \@off,
	  \@sdip, \@dep, \@idip, \@stop, \@reduction);
}

# last expression returns true value
1;
