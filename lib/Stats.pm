#!/usr/bin/perl -w

package Stats;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw( BINS MEAN STDEV ZSCORE LOG2);

sub BINS {
        my($binwidth,$bflank,$aflank)=@_;
        my @bincontainer;
        my $i;
        for ( $i=-$bflank;$i<=$aflank;$i+=$binwidth) {
                #next if $i == 0 ;
                push @bincontainer,$i ;
        }
        return @bincontainer ;
}

sub MEAN{
        my($data) = @_;
        if (not @$data) {
                die("Empty array $data\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = sprintf ("%.6f", $total / @$data ) ;
        return $average;
}
sub STDEV{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &MEAN($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = sprintf ( "%.6f", ($sqtotal / (@$data-1)) ** 0.5 ) ;
        my $err = sprintf ( "%.6f", ($std / sqrt(@$data)));
        return ($std,$err);
}

sub ZSCORE{
        my ($data,$stdev,$average)= @_;
        if ($stdev == 0 ){
                return (0);
        }
        my $zscore = ( $data - $average) / $stdev ;
        $zscore = sprintf ("%.6f", $zscore) ;
        return $zscore ;
}

sub ZSCORES{
	my($data)=@_;
	if(@$data == 1){
		return 0;
        }
	my @szcores;
	my $average = MEAN($data);
	my ($std,$err) = STDEV ($data);
	foreach (@$data){
		my $zscore=ZSCORE($_,$std,$average);
		push @szcores, $zscore ;
	}
	return(\@szcores);
}

sub SS {
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &MEAN($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        return ($sqtotal);
}

sub LOG2 {
       my $n = shift ;
	return log($n) / log(2);
}

1;
