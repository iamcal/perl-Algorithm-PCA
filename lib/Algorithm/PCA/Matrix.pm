package Algorithm::PCA::Matrix;

use base 'Math::MatrixReal';

use Data::Dumper;


sub new_ones {
	my $self = shift;
	$self = $self->new(@_);

	my ($m, $n) = $self->dim();

	for my $i (1..$m){
		for my $j (1..$n){
			$self->assign($i, $j, 1);
		}
	}

	return $self;
}

sub empirical_mean_rows {
	my ($self) = @_;

	my ($m, $n) = $self->dim();

	my $u = Algorithm::PCA::Matrix->new($m, 1);

	for my $i (1..$m){
		my $row = $self->row($i);
		my $sum = $row->norm_sum();
		$u->assign($i, 1, $sum / $n);
	}

	return $u;
}

sub ret_transpose {
	my ($self) = @_;

	my ($m, $n) = $self->dim();
	my $x =  Algorithm::PCA::Matrix->new($n, $m);
	$x->transpose($self);
	return $x;
}

sub covariance_dumb {
	my ($self) = @_;

	my ($m, $n) = $self->dim();

	my $b_trans = $self->ret_transpose();

	my $b_left = $self->shadow();
	$b_left->multiply_scalar($self, 1 / $n);

print $b_left."\n";
print $b_trans."\n";

	return $b_left * $b_trans;




	my $b_temp = $self * $b_trans;

	my $c = $b_temp->shadow();
	$c->multiply_scalar($b_temp, 1 / $m);

	return $c;
}


#
# i've confirmed that this one works correctly. woo
#

sub covariance {
	my ($self) = @_;

	# we'll assume that we have one dimension per row, so $n dims

	my ($n, $m) = $self->dim();


	#
	# get the means
	#

	my $means = [];
	for my $i (1..$n){
		$means->[$i] = 0;
		for my $j (1..$m){
			$means->[$i] += $self->element($i, $j);
		}
		$means->[$i] /= $m;
	}


	#
	# create each element of the covariance matrix
	#

	my $x = Algorithm::PCA::Matrix->new($n, $n);

	for my $i (1..$n){
		for my $j (1..$n){

			#
			# this is the covariance cov(i,j) of $self
			#

			my $total = 0;
			for my $k (1..$m){
				$total += ($self->element($i, $k) - $means->[$i]) * ($self->element($j, $k) - $means->[$j]);
			}

			$x->assign($i, $j, $total / ($m - 1));
		}
	}

	return $x;
}

sub clone_chunk {
	my ($self, $m, $n) = @_;

	my $x =  Algorithm::PCA::Matrix->new($m, $n);

	for my $i (1..$m){
		for my $j (1..$n){
			$x->assign($i, $j, $self->element($i, $j));
		}
	}

	return $x;
}

sub empirical_standard_deviation {
	my ($self) = @_;

	my ($m, $n) = $self->dim();

	my $x = Algorithm::PCA::Matrix->new($m , 1);
	for my $i (1..$m){
		$x->assign($i, 1, sqrt $self->element($i, $i));
	}

	return $x;
}

sub divide_by_element {
	my ($self, $x) = @_;

	return $self->each(sub {
		return 0 unless $x->element($_[1], $_[2]);
		$_[0] / $x->element($_[1], $_[2]);
	});
}

sub cumulative_energy_row {
	my ($self) = @_;

	my ($m, $n) = $self->dim();

	my $x = Algorithm::PCA::Matrix->new($m, 1);
	my $max = 0;

	for my $i (1..$m){
		my $total = 0;
		for my $j (1..$i){
			$total += $self->element($j, 1);
		}
		$x->assign($i, 1, $total);
		$max = $total;
	}

	for my $i (1..$m){
		$x->assign($i, 1, 100* $x->element($i, 1) / $max);
	}

	return $x;
}

1;
