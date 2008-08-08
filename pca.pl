#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;

my $x = [
	[1,0,0,1,0,1], # sample 1
	[0,2,0,1,0,2], # sample 2
	[1,0,1,0,1,0], # sample 3
];


#
# pre-step 1 - calculate matrix dims
#

my $n = scalar @{$x};
my $m = scalar @{$x->[0]};


#
# step 1 - calculate the empirical mean
#

my $u = [];

for my $i (1..$m){

	my $total = 0;

	for my $j (1..$n){
		$total += $x->[$j-1]->[$i-1];
	}

	push @{$u}, $total / $n;
}


#
# step 2 - calculate derivations from the mean
#

my $b = [];

for my $j (1..$n){

	my $row = [];

	for my $i (1..$m){
		push @{$row}, $x->[$j-1]->[$i-1] - $u->[$i-1];
	}

	push @{$b}, $row;
}


#
# step 3 - find the covariance matrix
#

my $b_div = matrix_multiply_scalar($b, 1 / $n);
my $b_trans = matrix_transpose($b);

#my $c = matrix_multiply($b_div, $b_trans);


#&dump_matrix($b);
&dump_matrix($b_div);
&dump_matrix($b_trans);



#print Dumper $b;
#&dump_matrix($x);
#&dump_matrix(matrix_transpose($x));




sub matrix_transpose {
	my ($x) = @_;
	my $y = [];

	my $n = scalar @{$x};
	my $m = scalar @{$x->[0]};

	for my $i (0..$m-1){
		$y->[$i] = [];
		for my $j (0..$n-1){
			push @{$y->[$i]}, $x->[$j]->[$i];
		}
	}

	return $y;
}

sub dump_matrix {
	my ($x) = @_;
	for my $row (@{$x}){
		for my $data (@{$row}){

			printf "% 02.2f ", $data;
		}
		print "\n";
	}
	print "\n";
}

sub matrix_multiply_scalar {
	my ($x, $s) = @_;

	my $n = scalar @{$x};
	my $m = scalar @{$x->[0]};
	my $y = [];

	for my $j (0..$n-1){
		$y->[$j] = [];
		for my $i (0..$m-1){
			$y->[$j]->[$i] = $x->[$j]->[$i] * $s;
		}
	}

	return $y;
}
