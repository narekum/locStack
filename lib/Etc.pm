package Etc;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw( STRANDCORRECTION ORIENTATION );




sub STRANDCORRECTION {
        my ($strand,$orientation)=@_;
        if ( $orientation eq "+") { return $strand ;}
        if ( $orientation eq "-" && $strand eq "fwd") { return "rev" ;}
        if ( $orientation eq "-" && $strand eq "rev") { return "fwd" ;}
        if ( $orientation eq "-" && $strand eq "both") { return "both" ;}
}


sub ORIENTATION {
        my ( $bin, $orientation )=@_;
        if ( $orientation eq "+" ) {
                return $bin ;
        } elsif ( $orientation eq "-" ) {
                return (-$bin) ;
        } else {
                return $bin ;
        }
}

