#!/usr/bin/perl -w

###############################################################################
# Splits list of coordinates from rhea into separate files, where each file
# contains coordinates at the same depth. Additionally, a list of all depth is
# written to a new file.
#
# Usage: perl split_coordinates.pl <input file>
###############################################################################

use strict;
use POSIX;

die "Usage: $0 <input file>\n" unless ($#ARGV == 0);

my $inputfile = $ARGV[0];

my $earth = 6371;

sub rheaz_to_depth {
    return $earth * (1. - $_[0]);
}

# The number of distinct z values encountered so far
my $no_files = 0;

# Store the file descriptor serial number for each encountered z value
my %z_hash;

# Store each file descriptor, depth and line count by their serial number
my @fhs;
my @fzs;
my @fls;

my $i = 0;

# Split input file and open outputs as necessary
open (IN, $inputfile) or die "Cannot open $inputfile\n";
while (<IN>) {
  my $origline = $_;
  chomp;
  #if (!($_ =~ m/(\d+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)/)) {
  if (!($_ =~ m/\s*(\d+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)\s+([-+\d\.e]+)/)) {
    die ("The following line does not match:\n$_\n");
  }
  my $rad0 = floor (1. * rheaz_to_depth ($2) + .5);
  my $rads = sprintf "%04d", $rad0;
  my ($outfile,$fh,$fhno);

  if (!exists $z_hash{$rad0}) {
    printf "Found new depth value no %d: %d\n", $no_files, $rad0;
    $fhno = $no_files;
    ($outfile = $inputfile) =~ s/\.txt$/_$rads.txt/;
    open ($fh, ">$outfile") or die "Cannot open $outfile\n";
    $z_hash{$rad0} = $fhno;
    $fhs[$fhno] = $fh;
    $fzs[$fhno] = $rad0;
    $fls[$fhno] = 1;
    ++$no_files;
  }
  else {
    $fhno = $z_hash{$rad0};
    $fh = $fhs[$fhno];
    ++$fls[$fhno];
  }
  print $fh $origline;
}

# Close all open files
printf "Found $no_files depths total\n";
for (my $i = 0; $i < $no_files; ++$i) {
  printf ("Found %d lines for depth %d\n", $fls[$i], $fzs[$i]);
  close $fhs[$i];
}
close($inputfile);

# Write depth list into file
my ($outfile,$fh,$rad0);

($outfile = $inputfile) =~ s/\.txt$/_z.txt/;
open ($fh, ">$outfile") or die "Cannot open $outfile\n";
foreach $rad0 (sort { $a > $b } @fzs) {
  printf $fh "%07d\n", $rad0;
}
close $fh;
