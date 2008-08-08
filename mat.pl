#!/usr/bin/perl

use warnings;
use strict;

use Math::MatrixReal;
use Data::Dumper;

#
# http://en.wikipedia.org/wiki/Principal_components_analysis
#

my $x = Math::MatrixReal->new_from_rows([
	[1,0,0,1,0,1], # sample 1
	[0,2,0,1,0,2], # sample 2
	[1,0,1,0,1,0], # sample 3
]);


#
# Calculate the empirical mean
#

my ($rows, $cols) = $x->dim();

my $u = Math::MatrixReal->new(1, $cols);

for my $i (1..$cols){
	my $col = $x->column($i);
	my $sum = $col->norm_sum();
	$u->assign(1, $i, $sum / $rows);
}


#
# Calculate the deviations from the mean
#

my $bx = $x->shadow();

for my $i (1..$rows){
	my $row = $x->row($i);
	my $new = $row - $u;
	$bx->assign_row($i, $new);
}


#
# rotate the dimensions - each sample is now in a column
#

$bx = &trans($bx);

my ($m, $n) = $bx->dim();

print "Dimensions (m) = $m\n";
print "Samples (n) = $n\n";

#&dump($bx);
#exit;


#
# Find the covariance matrix
#

my $b_trans = &trans($bx);

my $b_temp = $bx * $b_trans;

my $c = $b_temp->shadow();
$c->multiply_scalar($b_temp, 1 / $rows);

#&dump($c);
#exit;


#
# Find the eigenvectors and eigenvalues of the covariance matrix
#

my ($l, $v) = $c->sym_diagonalize();

my @pairs;

my $d = Math::MatrixReal->new($m, $m);
for my $i (1..$m){
	$d->assign($i, $i, $l->element($i, 1));
	push @pairs, [$i, $l->element($i, 1)];
}

#print Dumper $l->column(1);
#&dump($l);
#print Dumper \@pairs;
#&dump($d);
#&dump($v);
#exit;
#print "-------------------------------\n";


#
# Rearrange the eigenvectors and eigenvalues
#

my $d_prime = $d->shadow();
my $v_prime = $v->shadow();

my @sorted_pairs = sort { $b->[1] <=> $a->[1] } @pairs;

#print Dumper \@sorted_pairs;

my $col = 1;
for my $pair (@sorted_pairs){

	# set col $col in d_prime to col $pair->[0] in d
	# set col $col in v_prime to vol $pair->[0] in v

	for my $i (1..$m){

		$v_prime->assign($i, $col, $v->element($i, $pair->[0]));
	}

	$d_prime->assign($col, $col, $pair->[1]);

	$col++;
}

#&dump($d_prime);
#&dump($v_prime);
#exit;


#
# Compute the cumulative energy content for each eigenvector
#

my $g = Math::MatrixReal->new($m, 1);

for my $i (1..$m){
	my $total = 0;
	for my $j (1..$i){
		$total += $d_prime->element($j, $j);
	}
	$g->assign($i, 1, $total);
}

#&dump($g);


#
# Select a subset of the eigenvectors as basis vectors
#

my $limit = 2;

my $w = Math::MatrixReal->new($m, $limit);

for my $i (1..$m){
	for my $j (1..$limit){

		$w->assign($i, $j, $v->element($i, $j));
	}
}

#&dump($v);
#&dump($w);
#exit;


#
# Convert the source data to z-scores
#

my $s = Math::MatrixReal->new($m , 1);
for my $i (1..$m){
	$s->assign($i, 1, sqrt $c->element($i, $i));
}

my $h = Math::MatrixReal->new(1, $n);
for my $i (1..$n){
	$h->assign(1, $i, 1);
}

my $sh = $s * $h;

my $z = $bx->each(sub { $_[0] / $sh->element($_[1], $_[2]); });

#&dump($sh);
#&dump($bx);
#&dump($z);


#
# Project the z-scores of the data onto the new basis
#

my $w_trans = &trans($w);
my $y = $w_trans * $z;

&dump($y);




sub trans {
	my ($x) = @_;
	my ($rows, $cols) = $x->dim();
	my $new_x =  Math::MatrixReal->new($cols, $rows);
	$new_x->transpose($x);
	return $new_x;
}


sub dump {
	$_[0]->display_precision(2);
	print $_[0]."\n";
}
