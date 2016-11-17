##!/usr/bin/perl -w
# by Narendra Kumar @ 17/05/2016
#
# Apapted from Guangyao Li @ 03/10/2011
# 
#
###############################################################################################################
package GenomeInfo;
use strict;
use warnings;

use Exporter qw ( import );
our @EXPORT_OK = qw ( StdChromosomes );

sub StdChromosomes {
	my ($genome)=@_;
	my @chrom_information=();

	if($genome eq "dm2")
	{	push @chrom_information,	"chr2h", "chr2L", "chr2R", "chr3h", "chr3L", "chr3R", "chr4", "chr4h", "chrX", 
														"chrXh", "chrYh";
	}

	elsif($genome eq "dm3") 
	{	push @chrom_information,	"chr2L", "chr2LHet", "chr2R", "chr2RHet", "chr3L", "chr3LHet", "chr3R", "chr3RHet", 
														"chr4", "chrX", "chrXHet", "chrYHet";
	}

	elsif($genome eq "hg18")
	{	push @chrom_information,	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
														"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
														"chr20", "chr21", "chr22", "chrX", "chrY"; 
	}

	elsif($genome eq "hg19")
	{	push @chrom_information,	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
														"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
														"chr20", "chr21", "chr22", "chrX", "chrY"; 
	}

	elsif($genome eq "mm8")
	{	push @chrom_information,	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
														"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
														"chrX", "chrY"; 
	}

	elsif($genome eq "mm9")
	{	push @chrom_information,	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
														"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
														"chrX", "chrY"; 
	}

	elsif($genome eq "pombe")
	{	push @chrom_information,	"chr1", "chr2", "chr3", "mat"; 
	}

	elsif($genome eq "rn4")
	{	push @chrom_information,	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
														"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
														"chr20", "chrX"; 
	}

	elsif($genome eq "sacCer1")
	{	push @chrom_information,	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
														"chr11", "chr12", "chr13", "chr14", "chr15", "chr16"; 
	}

	else 
	{
		print " Error: Genome $genome not found\n";
	}
	return (\@chrom_information);
}
1;
