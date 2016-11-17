#!/usr/bin/perl

BEGIN {
        use Cwd qw(realpath cwd);
        use File::Basename;
        our ($fn, $dir) = fileparse(realpath($0));
}

use lib "$dir/../lib" ;

#use strict;
use LWP::Simple;
use Getopt::Long;
use Zdist;
use SeqToProfile;

my $cmd=$0." ".join(" ",@ARGV); ### command line copy

my $help;
my $matrixBed;
my $outfile;
my $binwidth;
my $bflank;
my $aflank;
my $scoretype="z";

GetOptions ('h|help'=>\$help,                     # --help         : print this help
            "m|mat=s" => \$matrixBed,             # --mat          : Cut matrix bed file produced by punarvasu.pl (required)
            "h|hex=s" => \$expmatrixBed,          # --hex          : Hexanucleotide matrix bed file (required)
#            "o|outfile=s" => \$outfile,          # --outfile      : outfile file name (required)
            "b|bin=i" => \$binwidth,              # --bin          : size of each bin in bp (required) 
            "d|down=i" => \$bflank,               # --down         : Window downstream of the reference coordinates (required)
            "u|up=i" => \$aflank,                 # --up           : Window upstream of the reference coordinates (required)
            "s|score=s" => \$scoretype,           # --score        : output mean(m) or z-scores(z) (default z)
            "p|preffile=s" => \$frequencyfile     # --preffile     : Nucleotide cut prefrences file
           )
or &PRINTHELP($fn);

if(defined $help || !$matrixBed || !$binwidth || !$bflank || !$aflank ) {
        &PRINTHELP($fn);
        exit;
}

####################################################################

my @options=( "$cmd",
           "--mat            $matrixBed",
#           "--outfile        $outfile",
           "--hex            $expmatrixBed",
           "--bin            $binwidth",
           "--down           $bflank",
           "--up             $aflank",
           "--scoretype      $scoretype",
           "--preffile       $frequencyfile"
);

print "** Running $fn with the following options\n\n" ;
print "      $_\n" foreach @options ;
print "\n";

# getting number of cuts 
my @matrixbed=Zdist::READMATRIXBED($matrixBed);
my ($matrix,$i,$j)=Zdist::CONVERTTOMATRIX(\@matrixbed);
my $stores=&Zdist::ZSCORESfromMATRIX($matrix,$i,$j);
my $totalofmeans=0;
foreach ( 1 .. $j){
	$totalofmeans+=$$stores[$_][1]; # [1] is mean; [4] is z-score
}
print "total of means=$totalofmeans\n";

#getting expected number of cuts
my @expmatrixbed=Zdist::READMATRIXBED($expmatrixBed);  # Hexanucleotide matrixbed around the center
my ($expmatrix,$i,$j)=Zdist::CONVERTTOMATRIX(\@expmatrixbed); #Hexanucleotide matrix around the center
my %freqhash=SeqToProfile::READFREQUENCIES($frequencyfile);
my $expectedcuts=Zdist::EXPECTEDCUTSPROFILE($expmatrix,$i,$j,\%freqhash);
#for (my $n=-$bflank;$n<=$aflank;$n+=$binwidth) {
#	my $i = (( $n - -$bflank ) / $binwidth) + 1  ;
#	my $exp_norm=$$expectedcuts[$i][2];
#	$totalexp_norm+=$exp_norm;
	#print "$n\t$exp_norm\t$totalexp_norm\n";
#}
my $scoreindex=1; # mean
#if ( $scoretype eq "Z" || $scoretype eq "z") { $scoreindex = 4 } 
#elsif ( $scoretype eq "M" || $scoretype eq "m") { $scoreindex = 1 }
#else { print "ERROR: unknown --score type\n" ; exit; };

$outfile=$matrixBed;
$outfile=~s/.mat//;
$outfile= $outfile . ".d" . $bflank . ".u" . $aflank . ".b" . $binwidth . ".s". $scoretype . ".ecr";
print "$outfile\n";
open ZWRITE,">$outfile";

my @expected_cuts_column=();
my @mean_cuts_column=();
for (my $n=-$bflank;$n<=$aflank;$n+=$binwidth) {
	my $i = (( $n - -$bflank ) / $binwidth) + 1  ;
	my $expectedcutscolumn=sprintf("%.6f", $$expectedcuts[$i][2] * $totalofmeans );
	push @expected_cuts_column,$expectedcutscolumn;
	push @mean_cuts_column, $$stores[$i][$scoreindex];
	#my $log2change=Stats::LOG2( $$stores[$i][$scoreindex] / $expectedcutscolumn ) if $expectedcutscolumn != 0 ;
	#print ZWRITE "$n", "\t",$$stores[$i][$scoreindex],"\t",$$expectedcuts[$i][1],"\t",$$expectedcuts[$i][2],"\t",$expectedcutscolumn ,"\t",$log2change,"\t",$$expectedcuts[$i][3], "\n";
}

$expected_cuts_column_zscore=Stats::ZSCORES(\@expected_cuts_column);
$mean_cuts_column_zscore=Stats::ZSCORES(\@mean_cuts_column);

for (my $n=-$bflank;$n<=$aflank;$n+=$binwidth) {
	my $i = (( $n - -$bflank ) / $binwidth) + 1  ;
	my $index=$i-1;
	my $log2change=Stats::LOG2( $$stores[$i][1] / $expected_cuts_column[$index] ) if $expected_cuts_column[$index] != 0 ;
	#my $log2change_zscore=Stats::LOG2( $$mean_cuts_column_zscore[$index] / $$expected_cuts_column_zscore[$index] ) if $$expected_cuts_column_zscore[$index] != 0 ;
	print ZWRITE "$n", "\t",$$stores[$i][1],"\t", $$mean_cuts_column_zscore[$index],"\t",$expected_cuts_column[$index],"\t",$$expected_cuts_column_zscore[$index],"\t",$log2change,"\t",$$expectedcuts[$i][2],"\n";
}


exit;

#for ( 1 .. 201 ) {
	#print "$stores[$_][1]\t$stores[$_][2]\t$stores[$_][3]\t$stores[$_][4]\n";
	#print "$$stores[$_][1]\t$$stores[$_][2]\t$$stores[$_][3]\t$$stores[$_][4]\t$$stores[$_][5]\t$$stores[$_][6]\n";
#	print "$$stores[$_][4]\n";
#}

#################################################################################################
sub PRINTHELP {
        my ($fn)=@_;

        print <<DES;

Usage: $fn < -m MATRIXBED -b BIN -d DOWNSTREM -u UPSTREAM > [Options]

Options:
  -h, --help       Show this help message.
  -m, --mat        Matrix file of reads in BED format as produced by punarvasu.pl program  (REQUIRED)
  -h  --hex        Hexanucleotide matrix bed file (required)
  -b, --bin        Size of each bin in bp 
  -d, --down       Window downstream of the reference coordinates. 
  -u, --up         Window upstream of the reference coordinates. 
  -p, --preffile   Prefrence file for the hexanucleotides

Examples: 
\$ perl $fn -mat MATRIXBED -bin BIN -down DOWNSTREM -up UPSTREAM

DES
exit;
}

