#!/usr/bin/perl -w

package ChattaIdx;
use strict;
use warnings;
use Fcntl;
use GenomeInfo;
use Parallel::ForkManager;
use CreateIndex;

use Exporter qw(import);

our @EXPORT_OK = qw( BINS SPTILREADHASH  READALIGNEDTAGS READSUMMITS STRANDCORRECTION ORIENTATION REPORTCOUNTSFWDSUMMIT REPORTCOUNTSREVSUMMIT);

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

sub READALIGNEDTAGS {
        my ($bedFile,$extend_length,$rscore,$dcol) =@_ ;
        my ($totalExpReads,@split,$chr,$each);
	my (@temp);
	my %alignedreadhash;
        open READS,$bedFile;
        print "      Loading tags from $bedFile File\n\n";
        #chomp (@temp=<READS>) ;
        #close READS ;
        #print STDERR "	READALIGNEDTAGS: Reading $bedFile File Read -- done\n";

        my $n=0;
        while ( $each = <READS> ) {
                next if $each =~ /^#/ ;
                chomp $each;
                @split=();
                @split=split(/\t/,$each);
                $chr=$split[0];
                $totalExpReads++;
                if ( $split[5] eq "-") {
                        if ( $rscore eq "n" ) {
                                $alignedreadhash{ $split[0] }{ $split[2] - 1 }{"-"}++;  # only records the last coord on the negative strand i.e the start position in 0 based.
                        } else {
                                $alignedreadhash{ $split[0] }{ $split[2] - 1 }{"-"}+=$split[$dcol - 1]  ;
                        }
                } else {
                        if ( $rscore eq "n" ){
                                $alignedreadhash{ $split[0] }{$split[1]}{"+"}++;
                        } else {
                                $alignedreadhash{ $split[0] }{$split[1]}{"+"}+=$split[$dcol - 1] ;
                        }
                }
                $n++;
                if ($n % 1000000 == 0 ){
			print "      No of tags processed...$n\n";
                }
        }
	print "\n";
	@temp=();

        return (\%alignedreadhash,$totalExpReads) ;
}

sub READALIGNEDTAGS_CHR {
        my ($bedFile,$extend_length,$rscore,$dcol,$chr_process) =@_ ;
	#print MINUS "$nuc\t$fractionminus\n";
        my ($totalExpReads,@split,$chr,$each);
        my (@temp);
        my %alignedreadhash;

	my $indexname = "$bedFile.index";
        my $indexoffset = "$bedFile.offset";

	my $chr_start_position, my $chr_end_position;
	open OPEN,$indexoffset or die "Can not open $indexoffset for read:$!\n";
	while (<OPEN>){
		chomp;
		my (@chroffset)=split(/\t/,$_);
		if ($chr_process eq $chroffset[0]){
			$chr_start_position=$chroffset[1];
			$chr_end_position=$chroffset[2];
			last;
		}
	}
	close OPEN ;
	print "$chr_process\t$chr_start_position\t$chr_end_position";
	#&CreateIndex::CREATE_OFFSET_INDEX($bedFile);
        #open READS,$bedFile;
        print "      Loading tags from $chr_process of $bedFile File\n\n";
        #chomp (@temp=<READS>) ;
        #close READS ;
        #print STDERR " READALIGNEDTAGS: Reading $bedFile File Read -- done\n";

        open(ORIG, "< $bedFile") or die "Can't open $bedFile for reading: $!";
        sysopen(IDX, $indexname, O_CREAT|O_RDWR,0666) or die "Can't open $indexname for read/write: $!";

        my $n=0;
        foreach ( $chr_start_position .. $chr_end_position ) {
		$each=CreateIndex::LINE_WITH_INDEX(*ORIG, *IDX,$_);
		#$each = CreateIndex::RETURN_LINE($bedFile,$_);
                chomp $each;
		#print "reading from index\t$each\t$n\n";
                next if $each =~ /^#/ ;
                @split=();
                @split=split(/\t/,$each);
                $chr=$split[0];
                $totalExpReads++;
                if ( $split[5] eq "-") {
                        if ( $rscore eq "n" ) {
                                $alignedreadhash{ $split[0] }{ $split[2] - 1 }{"-"}++;  # only records the last coord on the negative strand i.e the start position in 0 based.
                        } else {
                                $alignedreadhash{ $split[0] }{ $split[2] - 1}{"-"}+=$split[$dcol - 1]  ;
                        }
                } else {
                        if ( $rscore eq "n" ){
                                $alignedreadhash{ $split[0] }{$split[1]}{"+"}++;
                        } else {
                                $alignedreadhash{ $split[0] }{$split[1]}{"+"}+=$split[$dcol - 1] ;
                        }
                }
                $n++;
                if ($n % 1000000 == 0 ){
                        print "      No of tags processed on $chr_process...$n\n";
                }
        }
        print "\n";
        @temp=();

        return (\%alignedreadhash,$totalExpReads) ;
}



sub READSUMMITS {
        my ($summitFile)=@_;
        my (%summitsFwd,%summitsRev,$orientation,@temp,$chromosome,$start,$end,$totalsummits) ;
	my $totalfwdsummits=0;
	my $totalrevsummits=0;
	
        #push @groups,$group if $groupSize != 0 ;
        print "      Loading Reference coordinates from $summitFile File\n\n";
        open PEAKS,$summitFile;
        while (<PEAKS>){
                chomp;
                $totalsummits++;
                @temp=split(/\t/,$_ ) ;
                $start=$temp[1];
                $end=$temp[2];
                $chromosome=$temp[0];
		$orientation="+";
                if ( @temp == 3 ) {
                        $orientation="+";
		} elsif ($temp[5] eq "-" ) {
			$start = $end - 1 ;
			$orientation=$temp[5];
		} elsif ($temp[5] eq "+"){
			$orientation=$temp[5];
		}
		if ($orientation eq "+"){
			$totalfwdsummits++;
			$summitsFwd{$chromosome}{$start}=$orientation;
		} else {
			$totalrevsummits++;
			$summitsRev{$chromosome}{$start}=$orientation;
		}
        }
        close PEAKS ;
	print "      $totalsummits reference coordinates read from $summitFile\n";
	print "      $totalfwdsummits coordinates on Plus strand in $summitFile\n";
	print "      $totalrevsummits coordinates on Minus strand in $summitFile\n\n";
        return (\%summitsFwd,\%summitsRev,$totalfwdsummits,$totalrevsummits);
}
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
sub REPORTCOUNTSFWDSUMMIT {
        my ($alignedreads,$chr,$binstart,$binwidth,$extend_length) = @_ ;
        my $coordinates ;
        my $countboth=0 ;
	my $countminus=0;
	my $countplus=0;
        #if ($strand eq "fwd" ) { $strand = "+" } elsif ($strand eq "rev" ) { $strand = "-"; } #elsif { print "unknown strand type\n"; exit; }
        for $coordinates ( ( $binstart - $extend_length + 1 ) .. ( $binstart + $binwidth -1 )  ) {
                if (exists $$alignedreads{$chr}{$coordinates}{"+"} ) {
			#print "$chr $binstart +\n";
			$countplus += $$alignedreads{$chr}{$coordinates}{"+"} ;
                }
        }
	for $coordinates ( ( $binstart ) .. ( $binstart + $binwidth + $extend_length - 2 )  ) {
		if (exists $$alignedreads{$chr}{$coordinates}{"-"}  ){
			$countminus += $$alignedreads{$chr}{$coordinates}{"-"} ;
		}
	}
	$countboth=$countminus+$countplus;

    #    print "total counts= $count \n";
        return ($countplus,$countminus,$countboth) ;
}

sub REPORTCOUNTSREVSUMMIT {
        my ($alignedreads,$chr,$binstart,$binwidth,$extend_length) = @_ ;
        my $coordinates ;
        my $countboth=0 ;
        my $countminus=0;
        my $countplus=0;
        #if ($strand eq "fwd" ) { $strand = "+" } elsif ($strand eq "rev" ) { $strand = "-"; } #elsif { print "unknown strand type\n"; exit; }
        for $coordinates ( ( $binstart - $binwidth - $extend_length + 2 ) .. ( $binstart )  ) {
                if (exists $$alignedreads{$chr}{$coordinates}{"+"} ) {
                        #print "$chr $binstart +\n";
                        $countplus += $$alignedreads{$chr}{$coordinates}{"+"} ;
                }
        }
        for $coordinates ( ( $binstart - $binwidth +1 ) .. ( $binstart + $extend_length -  1 )  ) {
                if (exists $$alignedreads{$chr}{$coordinates}{"-"}  ){
                        $countminus += $$alignedreads{$chr}{$coordinates}{"-"} ;
                }
        }
        $countboth=$countplus+$countminus;

    #    print "total counts= $count \n";
        return ($countplus,$countminus,$countboth) ;
}
sub RETURNBINCOUNTS {
        my($summits,$alignedtagsexp,$alignedtagscontrol,$bincontainer,$binwidth,$extend_length,$normalize,$exptotalcounts,$controltotalcounts,$processchr,$rscore,$dcol)=@_;
	my($bin);
	my ($total_alignedtagsexp,$total_alignedtagscontrol);
	($alignedtagsexp,$total_alignedtagsexp)=READALIGNEDTAGS_CHR($alignedtagsexp,$extend_length,$rscore,$dcol,$processchr) if $alignedtagsexp ne "";
	($alignedtagscontrol,$total_alignedtagscontrol)=READALIGNEDTAGS_CHR($alignedtagscontrol,$extend_length,$rscore,$dcol,$processchr) if $alignedtagscontrol ne "NULL";

        #$normalize : per cent, per thousand, per million etc. The value should be 1000000 for per million
        #$exptotalcounts : total counts in exp
        #$controltotal counts : total couns in control exp
        my @bincontainer=@$bincontainer;
	my @returnchrmatrixBoth=();
	my @returnchrmatrixPlus=();
	my @returnchrmatrixMinus=();
	my @returnchrmatrix=();
        #print "--$_\t" foreach @bincontainer ; print "\n";
	foreach my $chr (keys %$summits) {
		next if $chr ne $processchr ;
		
		foreach my $summit ( keys $$summits{$chr}  ) {
			my $orientation=$$summits{$chr}{$summit};
			#print "$summit\n";
			#my @summitmatrix=(); 
			my (@summitmatrixPlus,@summitmatrixMinus,@summitmatrixBoth);
			my ($summitCorrectedCountBoth,$summitCorrectedCountPlus,$summitCorrectedCountMinus);
			push @summitmatrixBoth, $chr, $summit , $summit +1 , "both", "", $orientation ;
			push @summitmatrixPlus, $chr, $summit , $summit +1 , "plus", "", $orientation ;
			push @summitmatrixMinus, $chr, $summit , $summit +1 , "minus", "", $orientation ;
			foreach $bin ( @bincontainer ) {
				#print "t $bin "; print "\n";
				my ($countplusexp,$countminusexp,$countbothexp);
				my ($countpluscontrol,$countminuscontrol,$countbothcontrol);
				if ($orientation eq "-") {
					($countplusexp,$countminusexp,$countbothexp)=REPORTCOUNTSREVSUMMIT($alignedtagsexp,$chr,$summit - $bin,$binwidth,$extend_length);
				} else {
					($countplusexp,$countminusexp,$countbothexp)=REPORTCOUNTSFWDSUMMIT($alignedtagsexp,$chr,$summit + $bin,$binwidth,$extend_length);
				}
                                if ( $normalize != 0 ) {
                                        $countbothexp=sprintf ("%.1f", $countbothexp * $normalize / $exptotalcounts) ;
                                        $countplusexp=sprintf ("%.1f", $countplusexp * $normalize / $exptotalcounts) ;
                                        $countminusexp=sprintf ("%.1f", $countminusexp * $normalize / $exptotalcounts) ;
                                }

				if ($alignedtagscontrol ne "NULL"){
					if ($orientation eq "-") {
						($countpluscontrol,$countminuscontrol,$countbothcontrol)=REPORTCOUNTSREVSUMMIT($alignedtagscontrol,$chr,$summit - $bin,$binwidth,$extend_length);
					} else {
						($countpluscontrol,$countminuscontrol,$countbothcontrol)=REPORTCOUNTSFWDSUMMIT($alignedtagscontrol,$chr,$summit + $bin,$binwidth,$extend_length);
					}
	                                if ( $normalize != 0 ) {
                                	        $countbothcontrol=sprintf ("%.1f", $countbothcontrol * $normalize / $controltotalcounts) ;
                                        	$countpluscontrol=sprintf ("%.1f", $countpluscontrol * $normalize / $controltotalcounts) ;
                                        	$countminuscontrol=sprintf ("%.1f", $countminuscontrol * $normalize / $controltotalcounts) ;
                                	}
				} else { 
					$countpluscontrol=0,$countminuscontrol=0,$countbothcontrol=0;
				}

				my $correctedCountBoth=$countbothexp-$countbothcontrol;
				my $correctedCountPlus=$countplusexp-$countpluscontrol;
				my $correctedCountMinus=$countminusexp-$countminuscontrol;
				$summitCorrectedCountBoth+=$correctedCountBoth;
				$summitCorrectedCountPlus+=$correctedCountPlus;
				$summitCorrectedCountMinus+=$correctedCountMinus;	
				push @summitmatrixMinus,$correctedCountMinus;
				push @summitmatrixPlus,$correctedCountPlus;
				push @summitmatrixBoth,$correctedCountBoth;
				#print "$correctedCountBoth $correctedCountPlus $correctedCountMinus $chr $summit $orientation $bin\n";
			}
			#print "$_\t" foreach @summitmatrixPlus ; print "-\n";
			$summitmatrixPlus[4]=$summitCorrectedCountPlus;
			$summitmatrixMinus[4]=$summitCorrectedCountMinus;
			$summitmatrixBoth[4]=$summitCorrectedCountBoth;

			my $summitmatrixBoth_join=join("\t",@summitmatrixBoth); #print "-$summitmatrixBoth_join\n";
			my $summitmatrixPlus_join=join("\t",@summitmatrixPlus);
			my $summitmatrixMinus_join=join("\t",@summitmatrixMinus);

			push @returnchrmatrixBoth, join ("\t", @summitmatrixBoth);
			push @returnchrmatrixPlus, join ("\t", @summitmatrixPlus);
			push @returnchrmatrixMinus, join ("\t", @summitmatrixMinus);
			
		}
	}
	#print "$_\n" foreach @returnchrmatrixPlus ;
	@returnchrmatrix=(\@returnchrmatrixBoth,\@returnchrmatrixPlus,\@returnchrmatrixMinus);
	#return (\@returnchrmatrixBoth,\@returnchrmatrixPlus,\@returnchrmatrixMinus);
	return (\@returnchrmatrix);
}

sub RETURNBINCOUNTSALLCHR {
	my($summits,$alignedtagsexp,$alignedtagscontrol,$bincontainer,$binwidth,$extend_length,$normalize,$exptotalcounts,$controltotalcounts,$nproc,$genome,$rscore,$dcol)=@_;
	#my $chrmatrix;
	my @returnallchrmatrix=();
        my @returnallchrmatrixBoth=();
        my @returnallchrmatrixPlus=();
        my @returnallchrmatrixMinus=();
	my $chromosomes=GenomeInfo::StdChromosomes($genome);
	my @chromosomes=@$chromosomes;
	my (@bothtemp,@plustemp,@minustemp);
	#my $chrmatrix;

	my $pm = new Parallel::ForkManager($nproc);
 	$pm->run_on_finish(     # callback function to retrieve and assemble matrix from each chromosome matrices.
        sub{
                my($pid, $exit_code, $ident, $exit_signal, $core_dump, $dr) = @_;
                if(defined $dr){
                        print "      ## $ident \tFINISHED, distribution stored; pid: $pid\n";
                        my ($both,$plus,$minus) = @$dr ;
 		@bothtemp=();  @bothtemp=@$both;
 		@plustemp=();  @plustemp=@$plus;
 		@minustemp=(); @minustemp=@$minus;
 		push @returnallchrmatrixBoth,@bothtemp;
 		push @returnallchrmatrixPlus,@plustemp;
 		push @returnallchrmatrixMinus,@minustemp;
                        #@binmatrix{ keys %chrmatrix } = values %chrmatrix;
 		
                }else{
                        die "Child process returned void data in 'RETURNBINCOUNTSALLCHR' with exit code: $exit_code.\n";
                }
           }
        );

        $pm->run_on_start( sub {
                my ($pid,$ident)=@_;
                print "      ** $ident \tstarted, calculating distribution; pid: $pid\n";
        });
	
        foreach my $keys ( @chromosomes ) {
                my $pid = $pm->start($keys) and next;
                my $chrmatrix = ChattaIdx::RETURNBINCOUNTS ($summits,$alignedtagsexp,$alignedtagscontrol,$bincontainer,$binwidth,$extend_length,$normalize,$exptotalcounts,$controltotalcounts,$keys,$rscore,$dcol);
		#my $chrmatrix;
                $pm->finish(0, $chrmatrix);

        }

	$pm->wait_all_children;
	@returnallchrmatrix=(\@returnallchrmatrixBoth, \@returnallchrmatrixPlus, \@returnallchrmatrixMinus);
	print "\n";
	return (\@returnallchrmatrix);	
}
1;
