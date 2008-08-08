package Algorithm::PCA::Matrix;

use base 'Math::MatrixReal';



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

sub covariance {
	my ($self) = @_;

	my ($m, $n) = $self->dim();

	my $b_trans = $self->ret_transpose();

	my $b_temp = $self * $b_trans;

	my $c = $b_temp->shadow();
	$c->multiply_scalar($b_temp, 1 / $m);

	return $c;
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

	return $self->each(sub { $_[0] / $x->element($_[1], $_[2]); });
}

sub cumulative_energy_row {
	my ($self) = @_;

	 my ($m, $n) = $self->dim();

	my $x = Algorithm::PCA::Matrix->new($m, 1);

	for my $i (1..$m){
		my $total = 0;
		for my $j (1..$i){
			$total += $self->element($j, $j);
		}
		$x->assign($i, 1, $total);
	}
	return $x;
}

1;
