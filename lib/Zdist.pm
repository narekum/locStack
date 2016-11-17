#!/usr/bin/perl -w

package Zdist;
use strict;
use warnings;
use Stats;

use Exporter qw(import);

our @EXPORT_OK = qw( READMATRIXBED ZSCORESfromMATRIX PRINTMATRIX EXPECTEDCUTSPROFILE) ;

sub READMATRIXBED {
	my ($matrixBed) = @_ ;
	my @matrixbed=();
	open OPEN, $matrixBed or die "Could not open $matrixBed:$!\n";
	while (<OPEN>){
		next if $_=~/^\s*#/;
		chomp;
		push @matrixbed,$_;
	}
	return (@matrixbed) ;
}

sub CONVERTTOMATRIX {
	my ($matrixbed)=@_;
	my @matrix=();
	my $i=0; # row;
	my $j=0; # column;
	foreach my $row (@$matrixbed){
		$i++;
		my @split=split(/\t/,$row);
		for ( $j=6;$j<=$#split;$j++){
			$matrix[$i][$j-5]=$split[$j];
		}
	}
	return (\@matrix,$i,$j-6);
}

sub ZSCORESfromMATRIX {
	my ($matrix,$i,$j)=@_;
	my @allvalues=();
	my @storedscores=();
	my ($grandmean,$grandstdev,$grandsterr);
	#my $totalrows= length $#{$matrix};
	#print "total rws $i $j \n";
	foreach my $row ( 1 .. $i ) {
		foreach my $column ( 1 .. $j ) {
			#print "$$matrix[$row][$column] \n";
			push @allvalues,$$matrix[$row][$column];
		}
	}
	$grandmean=Stats::MEAN(\@allvalues);
	($grandstdev,$grandsterr)=Stats::STDEV (\@allvalues);
	#print "m s e -> $grandmean $grandstdev $grandsterr\n";

	foreach my $column ( 1 .. $j ) {
		my @binvalues=();
		foreach my $row ( 1 .. $i ) {
			push @binvalues, $$matrix[$row][$column] ;	
		}
		my $binmean=Stats::MEAN(\@binvalues);
		my ($binstdev,$binsterr)=Stats::STDEV (\@binvalues);
		#my ($binzscore)=Stats::ZSCORE($binmean,$binstdev,$grandmean);
		my $binzscore= sprintf ("%.4f", ( $binmean -$grandmean ) / $grandstdev );
		$storedscores[$column][0] = "" ; # in future may be add binname here. and column becomes row now.
		$storedscores[$column][1] = $binmean ; # 1 = mean ;
		$storedscores[$column][2] = $binstdev ; # 2 = stdev ;
		$storedscores[$column][3] = $binsterr ; # 3 = st err ;
		$storedscores[$column][4] = $binzscore ; # 4 = zscore ;
	}	
	return (\@storedscores);
}

sub EXPECTEDCUTSPROFILE {
	my ($matrix,$i,$j,$freqhash)=@_;
	my @storedscores=();
	my $sumOfExpectedScores=0;
	foreach my $column ( 1 .. $j ) {
		my %binvalues=();
		my $bintotal=0;
		my $expectedscore=0;
		foreach my $row ( 1 .. $i ) {
			$binvalues{$$matrix[$row][$column]}++; 
			$bintotal++;
		}
		foreach my $nuc ( keys %binvalues ) {
			$binvalues{$nuc} = sprintf ("%.6f",$binvalues{$nuc} / $bintotal) ;
		}
		#print "colu $column\n";
		#print "$_ \t$binvalues{$_}" , "\n" foreach keys %binvalues ;
		foreach ( keys %binvalues ) {
			$$freqhash{$_} = 0 if !$$freqhash{$_} ;
			$expectedscore += sprintf ("%.6f", $binvalues{$_} * $$freqhash{$_}) ; #print $binvalues{$_} , "-" , $$freqhash{$_} , "-" , $_ , "\n" if ! $$freqhash{$_} ; 
		}
		#print "colu $expectedscore\n";
		$storedscores[$column][1]=$expectedscore;
		$sumOfExpectedScores+=$expectedscore;
	}
	my @exp_freq_norm=();
	foreach my $column ( 1 .. $j ) {
		$storedscores[$column][2] = sprintf ("%.6f", $storedscores[$column][1] / $sumOfExpectedScores );
		push @exp_freq_norm,$storedscores[$column][2];
	}
	my $exp_freq_norm_mean=Stats::MEAN(\@exp_freq_norm);
	my $exp_freq_norm_std=Stats::STDEV(\@exp_freq_norm);
	foreach my $column ( 1 .. $j ) {
		my $exp_freq_norm_zscore=sprintf("%.6f", ( $storedscores[$column][2] - $exp_freq_norm_mean ) / $exp_freq_norm_std );
		$storedscores[$column][3]=$exp_freq_norm_zscore;
	}	
	
	return (\@storedscores);
	
}
sub PRINTMATRIX {
	my ($matrix,$i,$j,$outfile)=@_;
	open OPEN,">$outfile";
	foreach my $row ( 1 .. $i ) {
		foreach my $column ( 1 .. $j ) {
			print OPEN "$$matrix[$row][$column] \t";
			
		}
		print OPEN "\n";
	}
	close OPEN;
}

1;
