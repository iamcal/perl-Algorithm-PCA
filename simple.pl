#!/usr/bin/perl

use warnings;
use strict;

use lib 'lib';

use Algorithm::PCA;
use Algorithm::PCA::Matrix;
use Data::Dumper;

my $data = [

        [1,0,0,1,0,1], # sample 1
        [0,2,0,1,0,2], # sample 2
        [1,0,1,0,1,0], # sample 3

#[4,2,.6],
#[4.2,2.1,.59],
#[3.9,2.0,.58],
#[4.3,2.1,.62],
#[4.1,2.2,.63],

# http://www.maplepark.com/~drf5n/cgi-bin/c.cgi
#[1, 2, 3],
#[2, 4, 7],
#[1, 1, 2],
#[2, 3, 5],

#[3,0,1],
#[-4,1,2],
#[-6,0,-2],

#[2.5,0.5,2.2,1.9,3.1,2.3,2,1,1.5,1.1],     # i am a data series!
#[2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9],

];

my $out = Algorithm::PCA::reduce_simple($data);

my $comp = Algorithm::PCA::get_components($data);
my $out2 = Algorithm::PCA::reduce($data, $comp);

my $comp2 = Algorithm::PCA::get_components($data, 4);
my $out3 = Algorithm::PCA::reduce($data, $comp2);


print Dumper $out;
print Dumper $out2;
print Dumper $out3;
