#!/usr/bin/perl -w
use strict;
use POSIX;

# add `perl_libs` to the include path
use FindBin;
use lib "$FindBin::RealBin/perl_libs";

require Math::Spline;
require "myfunctions.pl";
require "readzones.pl";
use constant PI => 4*atan2(1,1);
use constant EPS => 1/10000;
use constant EARTH => 6371;
# turn autoflash on
$| = 1;

die "usage: ./weakzones inputfile depthfile between_depths outrefine\n"
  unless ($ARGV[0] and $ARGV[1]);
my $inputfile = $ARGV[0];
my $depthfile = $ARGV[1];


my $noskip = 1;
if ($ARGV[2]) {
  $noskip = $ARGV[2];
}

my $outrefine = 5;
if ($ARGV[3]) {
  $outrefine = $ARGV[3];
}

# read in depth file
my @depths;
my $no_depths;
open (DEP, $depthfile) or die "Cannot open $depthfile\n";
my $line = 0;
while (<DEP>) {
  if ($_ =~ m/\s?([\d\.e]+)\s?/) {
    # fill in intermediate depths
    my $nd = $1;
    if ($line > 0 and $depths[$line - 1] > $nd) {
      die "Depths in $depthfile need to be in ascending order\n";
    }
    if ($line > 0) {
      my $diff = $nd - $depths[$line - 1];
      for (my $k = 1; $k < $noskip; $k++) {
	$depths[$line] = $depths[$line - 1] + $diff / $noskip;
	$line++;
      }
    }
    $depths[$line] = $nd;
    $line++;
  } else {
    die "Cannot read line $line from $depthfile\n";
  }
}
$no_depths = $line;
close($depthfile);

my (@x1, @x2, @x3);
my (@x1n, @x2n, @x3n);
my ($k, $zone);
my (@param, @param_fine, @arc);
my $shift = 0;
my (@sdipn, @depn, @idipn, @offn, @stopn);
my $refine_factor = 1.1;
my $counter;

my ($nzones, $rnames, $rlnum, $rpvec, $rtvec, $roff,
    $rsdip, $rdep, $ridip, $rstop, $rreduction) = readtrench ($inputfile);

# de-reference and initialize
my @names = @$rnames;
my @lnum = @$rlnum;
my @pvec = @$rpvec;
my @tvec = @$rtvec;
my @off = @$roff;
my @sdip = @$rsdip;
my @dep = @$rdep;
my @idip = @$ridip;
my @stop = @$rstop;
my @reduction = @$rreduction;

# debug print
#  for ($zone = 0; $zone < $nzones; ++$zone) {
#       print "Data for ".$names[$zone]." \n";
#       my $npts = $lnum[$zone + 1] - $lnum[$zone];
#       my $offset = $lnum[$zone];
#       for (my $pts = 0; $pts < $npts; ++$pts) {
#   	printf "%3.2f %3.2f %.2f %.2f %.2f %.2f\n",
#   	    $pvec[$offset + $pts], $tvec[$offset + $pts],
#   	    $off[$offset + $pts], $sdip[$offset + $pts],
#   	    $dep[$offset + $pts], $idip[$offset + $pts];
#       }
#   }

my ($node0, $numnodes);

my $outfile_meta = sprintf("data/weakzones_meta.txt");
print "Opening  ".$outfile_meta.".\n";
open (OUTmeta, ">$outfile_meta") or die "Cannot open $outfile_meta\n";
print OUTmeta $nzones."\n";

for ($zone = 0; $zone < $nzones; ++$zone) {
  @x1 = []; @x2 = []; @x3 = [];
  @x1n = []; @x2n = []; @x3n = [];
  @arc = []; @param  = []; @param_fine = [];

  # cut points for one zone off
  $node0 = splice(@lnum,0,1);
  $numnodes = $lnum[0] - $node0;
  print "Number of nodes: ".$numnodes."\n";
  my @p = splice(@pvec,0,$numnodes);
  my @t = splice(@tvec,0,$numnodes);
  my @o  = splice(@off,0,$numnodes);
  my @s  = splice(@sdip,0,$numnodes);
  my @d  = splice(@dep,0,$numnodes);
  my @i  = splice(@idip,0,$numnodes);
  my @st  = splice(@stop,0,$numnodes);

  my $num_pt = floor($refine_factor * $numnodes);
  my $depth = 0;
  my $depth_old = 0;
  my ($dist, $t);
  my ($outfile, $outfile_xyz);
  my ($spx1, $spx2, $spx3, $spsdip, $spdep, $spidip, $spoff, $spstop);
  my ($xmin, $xmax, $ymin, $ymax, $zmin, $zmax);

  # convert into Carthesian coordinates on r = 1
  # and derive parametrization
  for ($k=0; $k<$numnodes; ++$k) {
    ($x1[$k], $x2[$k], $x3[$k]) = S2C (180 + $p[$k],90 - $t[$k]-$shift);
  }
  $arc[0] = 0;
  for ($k=1; $k<$numnodes; ++$k) {
    $dist = vec_norm ($x1[$k] - $x1[$k-1],
		      $x2[$k] - $x2[$k-1],
		      $x3[$k] - $x3[$k-1]);
    $arc[$k] = $arc[$k-1] + $dist;
  }
  # shift parameters to interval [0,1]:
  for ($k=0; $k<$numnodes; ++$k) {
    $arc[$k] /= $arc[$numnodes-1];
  }
  for ($k=0; $k<$num_pt; ++$k) {
    $param[$k] = $k/($num_pt - 1);
  }
  for ($k=0; $k<$outrefine*$num_pt; ++$k) {
    $param_fine[$k] = $k/($outrefine*$num_pt - 1);
  }

  # interplate data to finer mesh
  $spx1 = new Math::Spline(\@arc,\@x1);
  $spx2 = new Math::Spline(\@arc,\@x2);
  $spx3 = new Math::Spline(\@arc,\@x3);
  $spsdip = new Math::Spline(\@arc,\@s);
  $spidip = new Math::Spline(\@arc,\@i);
  $spdep = new Math::Spline(\@arc,\@d);
  $spoff = new Math::Spline(\@arc,\@o);
  $spstop = new Math::Spline(\@arc,\@st);

  # reparametrize to uniform t-values
  for ($k=0; $k<$num_pt; ++$k) {
    $t = $param[$k];
    $x1n[$k] = $spx1->evaluate($t);
    $x2n[$k] = $spx2->evaluate($t);
    $x3n[$k] = $spx3->evaluate($t);
    $sdipn[$k] = $spsdip->evaluate($t);
    $idipn[$k ] = $spidip->evaluate($t);
    $depn[$k] = $spdep->evaluate($t);
    $offn[$k] = $spoff->evaluate($t);
    $stopn[$k] = $spstop->evaluate($t);
  }
  @x1 = @x1n; @x2 = @x2n; @x3 = @x3n;
  @s = @sdipn; @i = @idipn; @d = @depn; @o = @offn; @st = @stopn;

  my ($ll, $x1f, $x2f, $x3f, $stopf, $p1, $t1);
  $outfile_xyz = sprintf("data/".$names[$zone]."_".$zone."_".$num_pt."_weak.xyz");
  print "Opening outfile ".$outfile_xyz.".\n";
  open (OUTxyz, ">$outfile_xyz") or die "Cannot open $outfile_xyz\n";
  $counter = 0;
  # loop over depths
  for ($ll = 0; $ll <= $no_depths; ++$ll) {
    if ($ll < $no_depths) {
      $depth = $depths[$ll]/EARTH;
    }

    # choose approximative arc length as parametrization
    $arc[0] = 0;
    for ($k=1; $k<$num_pt; ++$k) {
      $dist = vec_norm ($x1[$k] - $x1[$k-1],
			$x2[$k] - $x2[$k-1],
			$x3[$k] - $x3[$k-1]);
      $arc[$k] = $arc[$k-1] + $dist;
    }

    # shift parameters to interval [0,1]:
    for ($k=0; $k<$num_pt; ++$k) {
      $arc[$k] /= $arc[$num_pt-1];
    }

    # define splines
    $spx1 = new Math::Spline(\@arc,\@x1);
    $spx2 = new Math::Spline(\@arc,\@x2);
    $spx3 = new Math::Spline(\@arc,\@x3);
    $spsdip = new Math::Spline(\@arc,\@s);
    $spidip = new Math::Spline(\@arc,\@i);
    $spdep = new Math::Spline(\@arc,\@d);
    $spoff = new Math::Spline(\@arc,\@o);
    $spstop = new Math::Spline(\@arc,\@st);

    # reparametrize to uniform t-values
    for ($k=0; $k<$num_pt; ++$k) {
      $t = $param[$k];
      $x1n[$k] = $spx1->evaluate($t);
      $x2n[$k] = $spx2->evaluate($t);
      $x3n[$k] = $spx3->evaluate($t);
      $sdipn[$k] = $spsdip->evaluate($t);
      $depn[$k] = $spdep->evaluate($t);
      $idipn[$k] = $spidip->evaluate($t);
      $offn[$k] = $spoff->evaluate($t);
      $stopn[$k] = $spstop->evaluate($t);
    }

    # write slab in a file
    if ($ll == 0) {
      $outfile = sprintf("data/".$names[$zone]."_".$zone."_".$num_pt."_weak_%04.0f_data.txt", $depth*EARTH);
    } else {
      $outfile = sprintf("data/".$names[$zone]."_".$zone."_".$num_pt."_weak_%04.0f.txt", $depth_old*EARTH);
    }
    if ($ll == 0 || ((($ll - 1) % $noskip) == 0)) {
      print "Writing to outfile ".$outfile.".\n";
      open (OUT, ">$outfile") or die "Cannot open $outfile\n";
      # interplate to higher reslution for output
      for ($k=0; $k<$outrefine*$num_pt; ++$k) {
	$t = $param_fine[$k];
	$x1f = $spx1->evaluate($t);
	$x2f = $spx2->evaluate($t);
	$x3f = $spx3->evaluate($t);
	$stopf = $spstop->evaluate($t);
	if ($depth * EARTH <= $stopf) {
	  ($p1,$t1) = C2S ($x1f, $x2f, $x3f);
	  $t1 = 90 - $t1;
	  $p1 = $p1 - 180;
	  print OUT $p1+$shift."   ".$t1."\n";
	  if ($ll > 0) {
	    print OUTxyz $x1f."   ".$x2f."   ".$x3f." \n";
	    $counter++;
	  }
	  if ($k == 0 && $ll == 0) {
	    $xmin = $xmax = $x1f;
	    $ymin = $ymax = $x2f;
	    $zmin = $zmax = $x3f;
	  } else {
	    $xmin = min($xmin, $x1f);
	    $xmax = max($xmax, $x1f);
	    $ymin = min($ymin, $x2f);
	    $ymax = max($ymax, $x2f);
	    $zmin = min($zmin, $x3f);
	    $zmax = max($zmax, $x3f);
	  }
	}
      }
      close OUT;
    }

    # derive normal vector for each point
    for ($k=0; $k<$num_pt; ++$k) {
      # derive curve tangential
      my ($dx1, $dx2, $dx3) = curve_derivative ($spx1, $spx2, $spx3, $param[$k]);

      # derive unit length normal vector
      my ($n1,$n2,$n3) = cross_prod ($x1n[$k],$x2n[$k],$x3n[$k],
				     $dx1,$dx2,$dx3);
      my $nnorm = vec_norm ($n1,$n2,$n3);
      $n1 /= $nnorm;
      $n2 /= $nnorm;
      $n3 /= $nnorm;
      # surface shift by off in normal direciton
      if ($ll == 0) {
        $x1n[$k] += $offn[$k]/EARTH * $n1;
	$x2n[$k] += $offn[$k]/EARTH * $n2;
	$x3n[$k] += $offn[$k]/EARTH * $n3;
      } else {
	my $d1old = weakzone_fun($depth_old, $depn[$k]/EARTH,$sdipn[$k], $idipn[$k]);
	my $d1 = weakzone_fun($depth, $depn[$k]/EARTH, $sdipn[$k], $idipn[$k]);
	$x1n[$k] = ($d1 - $d1old) * $n1 + (1 + ($depth_old - $depth)) * $x1n[$k];
	$x2n[$k] = ($d1 - $d1old) * $n2 + (1 + ($depth_old - $depth)) * $x2n[$k];
	$x3n[$k] = ($d1 - $d1old) * $n3 + (1 + ($depth_old - $depth)) * $x3n[$k];
      }
    }
    @x1 = @x1n; @x2 = @x2n; @x3 = @x3n;
    @s = @sdipn; @d = @depn; @i = @idipn; @st = @stopn;
    $depth_old = $depth;
  }				# done with loop over depths
  close OUTxyz;
  print OUTmeta $names[$zone]."_".$zone."_".$num_pt."_weak.xyz\n";
  print OUTmeta $counter."\n";
  print OUTmeta $reduction[$zone]."\n";
  # put 10km buffer around point data
  $xmin -= 10.0/EARTH;
  $xmax += 10.0/EARTH;
  $ymin -= 10.0/EARTH;
  $ymax += 10.0/EARTH;
  $zmin -= 10.0/EARTH;
  $zmax += 10.0/EARTH;
  print OUTmeta $xmin."\n".$xmax."\n".$ymin."\n".
    $ymax."\n".$zmin."\n".$zmax."\n";
}
close OUTmeta;


