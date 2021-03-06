#!/usr/bin/perl -w

package SeqToProfile;
#use lib "/home/naren/build/punarvasu/lib" ;
use strict;
use warnings;
use GenomeInfo;

use Exporter qw(import);

our @EXPORT_OK = qw( READFREQUENCIES );

sub READFREQUENCIES {
	my ($freqfile)=@_ ;
	open OPEN,$freqfile or die ("Can not open $freqfile:$!\n");
	my %freqhash;
	while (<OPEN>){
		chomp;
		my @split=split(/\t/,$_);
		$freqhash{$split[0]}=$split[1];
	}
	return (%freqhash);
}

sub SeqToPref {
	my ($seq,$model)=@_; # model=6(hexamer),4(tetramer)
	my $length=length($seq);
	#print "seqlength-$seq\t$length\n";
	my @prefbed=();
	my @prefbed_oppo=();
	for (my $i=0;$i<=$length-$model;$i++) {
		my $nmer=substr($seq,$i,$model);
		my $nmerrevcom=REVCOM($nmer);
		#my $nmaervalue=$$freqhash{$nmer};
		#print "$nmer-$nmaervalue\t";
		push @prefbed,$nmer;
		push @prefbed_oppo,$nmerrevcom;
	}
	#print "\n";
	return (\@prefbed,\@prefbed_oppo);
}

sub BedToPref {
	my ($summits,$bflank,$aflank,$model,$genome,$fasta)=@_;
	#my $chromosomes=GenomeInfo::StdChromosomes($genome);
	open BEDWRITE,">temp.bed" or die "Cant not open temp for writing:$!\n";
	foreach my $chr (keys %$summits) {
		foreach my $summit ( sort {$a<=>$b} keys $$summits{$chr}  ) {
			my $orientation=$$summits{$chr}{$summit};
			my ($startentended,$endextended);
			if ($orientation eq "+") {
				$startentended = $summit - $bflank -  $model / 2 ;
				$endextended = $summit + $aflank + $model / 2  ;
			} elsif ( $orientation eq "-" ) {
				$endextended = $summit + $bflank + $model / 2 + 1  ;
				$startentended = $summit - $aflank -  $model / 2 + 1  ;
			}
			my $summitheader = join "\t" , $chr , $startentended , $endextended , "." , "0", $orientation ; 
			print BEDWRITE "$summitheader\n";
			#print "$summitheader\n";
		}
	}
	#close BEDWRITE ;
	system ( "bedtools getfasta -fi $fasta -bed temp.bed -fo temp.seq -tab -s" );
	open OPEN, "temp.seq";
	my @matrixout=();
	my @matrixout_oppo=();
	while (<OPEN>){
		chomp;
		my ($chr,$start,$end,$strand,$nuc)=/^(\S+):(\d+)-(\d+)\((\S)\)\t(\S+)$/;
		my ($prefbed, $prefbed_oppo )= SeqToPref ($nuc, $model);
		if ($strand eq "+") {
                	$start = $start + $bflank +  $model / 2 ;
                	$end = $start + 1 ;
                } elsif ( $strand eq "-" ) {
                        $end = $end - $bflank - $model / 2 - 1  ;
                        $start = $end - 1 ;
               	}
		my $oppo_first = shift @$prefbed_oppo ;
		my $oppo_last = pop @$prefbed_oppo ;
		#$oppo_last = substr ($oppo_last,0,5) . "A" ;
		push @$prefbed_oppo, $oppo_last , $oppo_last ;

		my $summitPrefs=join "\t" , $chr,$start,$end,"same","0",$strand, @$prefbed ;
		my $summitPrefs_oppo=join "\t" , $chr,$start,$end,"oppo","0",$strand, @$prefbed_oppo ;
		#print "$summitPrefs_oppo\n";
		push @matrixout, $summitPrefs ;
		push @matrixout_oppo, $summitPrefs_oppo ;
	}
	system ("rm temp.bed temp.seq");
	return (\@matrixout,\@matrixout_oppo);
}

sub REVCOM {
        my ($motif)=@_;
        $motif =~ tr/ACGTNacgtn/TGCANTGCAN/ ;
        $motif = reverse $motif ;
        return ($motif);
}

1;
