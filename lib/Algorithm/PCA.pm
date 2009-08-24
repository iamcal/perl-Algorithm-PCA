package Algorithm::PCA;

use warnings;
use strict;

our $VERSION = 0.01;

use Algorithm::PCA::Matrix;


#
# returns a components matrix which can be used to do PCA on a dataset
# (requires a dataset to derive this from)
#

sub get_components {
	my ($data, $dims) = @_;

	my $x = Algorithm::PCA::Matrix->new_from_cols($data);


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

	bless $d_prime, 'Algorithm::PCA::Matrix';
	bless $v_prime, 'Algorithm::PCA::Matrix';

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


	#&dump($d_prime);
	#&dump($v_prime);
	#exit;


	#
	# we can manually specify the dims, else figure it out based on energy
	#

	my $limit = $dims;

	if (!$limit){

		#
		# Compute the cumulative energy content for each eigenvector
		#

		my $g = $d_prime->cumulative_energy_row();


		#
		# we want enough dims to include 90%+ of the eigenvalues
		#

		for my $i (1..$m){
			$limit = $i;
			last unless $g->element($i, 1) < 90;
		}
	}

	#print "limit is $limit\n";
	#&dump($g);
	#exit;


	#
	# Select a subset of the eigenvectors as basis vectors
	#

	my $w = $v_prime->clone_chunk($m, $limit);


	return $w->ret_transpose();
}


#
# reduce a data set
#

sub reduce {
	my ($data, $components) = @_;

	my $x = Algorithm::PCA::Matrix->new_from_cols($data);

	my $reduced = $components * $x;


	#
	# convert into an arrayref of arrayrefs;
	#

	my $out = [];

	$reduced->each( sub { $out->[$_[2]-1]->[$_[1]-1] = $_[0]; } );

	return $out;
}


#
# reduce simple!
#

sub reduce_simple {
	my ($data) = @_;

	my $components = &get_components($data);
	my $out = &reduce($data, $components);

	return $out;
}






sub bleh {

	#
	# Convert the source data to z-scores
	#

#	my $z = $bx->divide_by_element($c->empirical_standard_deviation() * $h);

	#&dump($z);


	#
	# Project the z-scores of the data onto the new basis
	#

#	my $y = $w->ret_transpose() * $z;

#	return $y;
}



sub dump {
	$_[0]->display_precision(2);
	print $_[0]."\n";
}

1;
