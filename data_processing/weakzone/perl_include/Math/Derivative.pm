# functions for calculating derivatives of data
# $Id: Derivative.pm,v 1.1 1995/12/26 16:26:59 willijar Exp $
=head1 NAME

 Math::Derivative - Numeric 1st and 2nd order differentiation

=head1 SYNOPSIS

    use Math::Derivative qw(Derivative1 Derivative2);
    @dydx=Derivative1(\@x,\@y);
    @d2ydx2=Derivative2(\@x,\@y);
    @d2ydx2=Derivative2(\@x,\@y,$yp0,$ypn);

=head1 DESCRIPTION

This Perl package exports functions for performing numerical first
(B<Derivative1>) and second B<Derivative2>) order differentiation on
vectors of data. They both take references to two arrays containing
the x and y ordinates of the data and return an array of the 1st or
2nd derivative at the given x ordinates. B<Derivative2> may optionally
be given values to use for the first dervivative at the start and end
points of the data - otherwiswe 'natural' values are used.

=head1 HISTORY

$Log: Derivative.pm,v $
Revision 1.1  1995/12/26 16:26:59  willijar
Initial revision


=head1 BUGS

Let me know.

=head1 AUTHOR

John A.R. Williams <J.A.R.Williams@aston.ac.uk>

=cut

require Exporter;
package Math::Derivative;
@ISA=qw(Exporter);
@EXPORT_OK=qw(Derivative1 Derivative2);
use strict;

sub Derivative1 {
    my ($x,$y)=@_;
    my @y2;
    my $n=$#{$x};
    $y2[0]=($y->[1]-$y->[0])/($x->[1]-$x->[0]);
    $y2[$n]=($y->[$n]-$y->[$n-1])/($x->[$n]-$x->[$n-1]);
    my $i;
    for($i=1; $i<$n; $i++) {
	$y2[$i]=($y->[$i+1]-$y->[$i-1])/($x->[$i+1]-$x->[$i-1]);
    }
    return @y2;
}

sub Derivative2 {
    my ($x,$y,$yp1,$ypn)=@_;
    my $n=$#{$x};
    my ($i,@y2,@u);
    if (!defined $yp1) {
	$y2[0]=0; $u[0]=0;
    }
    else {
	$y2[0]=-0.5;
	$u[0]=(3/($x->[1]-$x->[0]))*(($y->[1]-$y->[0])/($x->[1]-$x->[0])-$yp1);
    }
    for($i=1; $i<$n; $i++) {
	my $sig=($x->[$i]-$x->[$i-1])/($x->[$i+1]-$x->[$i-1]);
	my $p=$sig*$y2[$i-1]+2.0; 	
	$y2[$i]=($sig-1.0)/$p;
	$u[$i]=(6.0*( ($y->[$i+1]-$y->[$i])/($x->[$i+1]-$x->[$i])-
		      ($y->[$i]-$y->[$i-1])/($x->[$i]-$x->[$i-1])
		     )/
		($x->[$i+1]-$x->[$i-1])-$sig*$u[$i-1])/$p;
    }
    my ($qn,$un);
    if (!defined $ypn) {
	$qn=0;
	$un=0;
    }
    else {
	$qn=0.5;
	$un=(3.0/($x->[$n]-$x->[$n-1]))*
	    ($ypn-($y->[$n]-$y->[$n-1])/($x->[$n]-$x->[$n-1]));
    }
    $y2[$n]=($un-$qn*$u[$n-1])/($qn*$y2[$n-1]+1.0);
    for($i=$n-1; $i>=0; --$i) {
	$y2[$i]=$y2[$i]*$y2[$i+1]+$u[$i];
    }
    return @y2;
}


1;
