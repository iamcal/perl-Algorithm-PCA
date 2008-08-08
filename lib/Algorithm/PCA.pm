package Algorithm::PCA;

use warnings;
use strict;

use Algorithm::PCA::Matrix;
use Data::Dumper;

sub reduce {
	my ($x) = @_;

	#
	# http://en.wikipedia.org/wiki/Principal_components_analysis
	#

	my ($m, $n) = $x->dim();

	#print "Dimensions (m) = $m\n";
	#print "Samples (n) = $n\n";
	#&dump($x);


	#
	# Calculate the empirical mean
	#

	my $u = $x->empirical_mean_rows();


	#
	# Calculate the deviations from the mean
	#

	my $h = Algorithm::PCA::Matrix->new_ones(1, $n);

	my $bx = $x - ($u * $h);


	#
	# Find the covariance matrix
	#

	my $c = $bx->covariance();


	#
	# Find the eigenvectors and eigenvalues of the covariance matrix
	#

	my ($d, $v) = $c->sym_diagonalize();


	#
	# @pairs contains [column_num, eigenvalue] pairs for sorting
	#

	my @pairs = map { [$_, $d->element($_, 1)] }(1..$m);


	#
	# Rearrange the eigenvectors and eigenvalues
	#

	my $d_prime = $d->shadow();
	my $v_prime = $v->shadow();

	my @sorted_pairs = sort { $b->[1] <=> $a->[1] } @pairs;

	my $col = 1;
	for my $pair (@sorted_pairs){

		# set col $col in d_prime to col $pair->[0] in d
		# set col $col in v_prime to vol $pair->[0] in v

		for my $i (1..$m){

			$v_prime->assign($i, $col, $v->element($i, $pair->[0]));
		}

		$d_prime->assign($col, 1, $pair->[1]);

		$col++;
	}


	#
	# Compute the cumulative energy content for each eigenvector
	#

	#my $g = $d_prime->cumulative_energy_row();

	#&dump($g);


	#
	# Select a subset of the eigenvectors as basis vectors
	#

	my $limit = 2;
	my $w = $v_prime->clone_chunk($m, $limit);


	#
	# Convert the source data to z-scores
	#

	my $z = $bx->divide_by_element($c->empirical_standard_deviation() * $h);

	#&dump($z);


	#
	# Project the z-scores of the data onto the new basis
	#

	my $y = $w->ret_transpose() * $z;

	return $y;
}



sub dump {
	$_[0]->display_precision(2);
	print $_[0]."\n";
}

1;
