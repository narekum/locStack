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
use Stats;

my $cmd=$0." ".join(" ",@ARGV); ### command line copy

my $help;
my $matrixBed;
my $expmatrixBed;
my $outfile;
my $binwidth;
my $bflank;
my $aflank;
#my $scoretype="z";

GetOptions ('h|help'=>\$help,                     # --help         : print this help
            "m|mat=s" => \$matrixBed,             # --mat          : Cut matrix bed file produced by punarvasu.pl (required)
            "e|expected=s" => \$expmatrixBed,        # --expected     : Expexted cut matrix bed file produced by punarvasu.pl (required)
#            "o|outfile=s" => \$outfile,           # --outfile      : outfile file name (required)
            "b|bin=i" => \$binwidth,              # --bin          : size of each bin in bp (required) 
            "d|down=i" => \$bflank,               # --down         : Window downstream of the reference coordinates (required)
            "u|up=i" => \$aflank,                 # --up           : Window upstream of the reference coordinates (required)
           )
or &PRINTHELP($fn);

if(defined $help || !$matrixBed || !$expmatrixBed || !$binwidth || !$bflank || !$aflank ) {
        &PRINTHELP($fn);
        exit;
}

####################################################################

my @options=( "$cmd",
           "--mat            $matrixBed",
           "--expected       $expmatrixBed",
#           "--outfile        $outfile",
           "--bin            $binwidth",
           "--down           $bflank",
           "--up             $aflank",
);

print "** Running $fn with the following options\n\n" ;
print "      $_\n" foreach @options ;
print "\n";

my @matrixbed=Zdist::READMATRIXBED($matrixBed);
my ($matrix,$i,$j)=Zdist::CONVERTTOMATRIX(\@matrixbed);
my $stores=&Zdist::ZSCORESfromMATRIX($matrix,$i,$j);

my @expmatrixbed=Zdist::READMATRIXBED($expmatrixBed);
my ($expmatrix,$i,$j)=Zdist::CONVERTTOMATRIX(\@expmatrixbed);
my $expstores=&Zdist::ZSCORESfromMATRIX($expmatrix,$i,$j);

#my $scoreindex;
#if ( $scoretype eq "Z" || $scoretype eq "z") { $scoreindex = 4 } 
#elsif ( $scoretype eq "M" || $scoretype eq "m") { $scoreindex = 1 }
#else { print "ERROR: unknown --score type\n" ; exit; };

$outfile=$matrixBed;
$outfile=~s/.mat//;
$outfile= $outfile . ".d" . $bflank . ".u" . $aflank . ".b" . $binwidth . ".log2.txt";
print "$outfile\n";
open ZWRITE,">$outfile";

for (my $n=-$bflank;$n<=$aflank;$n+=$binwidth) {
	my $i = (( $n - -$bflank ) / $binwidth) + 1  ;
        my $log2change = Stats::LOG2( $$stores[$i][1] / $$expstores[$i][1]) if $$expstores[$i][1] != 0 ;
	$log2change=5 if $$expstores[$i][1] == 0 ;
	print ZWRITE "$n", "\t",$log2change,"\n";
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
  -e, --expcutmat  Expected cuts Matrix file of reads in BED format as produced by punarvasu.pl program  (REQUIRED)
  -b, --bin        Size of each bin in bp 
  -d, --down       Window downstream of the reference coordinates. 
  -u, --up         Window upstream of the reference coordinates. 

Examples: 
\$ perl $fn -mat MATRIXBED -expcutmat EXPECTEDCUTSMATRIX -bin BIN -down DOWNSTREM -up UPSTREAM

DES
exit;
}

