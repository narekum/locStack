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
#            "o|outfile=s" => \$outfile,           # --outfile      : outfile file name (required)
            "b|bin=i" => \$binwidth,              # --bin          : size of each bin in bp (required) 
            "d|down=i" => \$bflank,               # --down         : Window downstream of the reference coordinates (required)
            "u|up=i" => \$aflank,                 # --up           : Window upstream of the reference coordinates (required)
            "s|score=s" => \$scoretype,           # --score        : output mean(m) or z-scores(z) (default z)
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
           "--bin            $binwidth",
           "--down           $bflank",
           "--up             $aflank",
           "--scoretype      $scoretype",
);

print "** Running $fn with the following options\n\n" ;
print "      $_\n" foreach @options ;
print "\n";

my @matrixbed=Zdist::READMATRIXBED($matrixBed);
#print "$_\n" foreach @matrixbed ;

$outfile=$matrixBed;
$outfile=~s/.mat/.mat.dat/;

my ($matrix,$i,$j)=Zdist::CONVERTTOMATRIX(\@matrixbed);
#print "$$matrix[1][$_]\t" foreach ( 1 .. 401); print "\n";

&Zdist::PRINTMATRIX($matrix,$i,$j,$outfile);

#my $stores=&Zdist::ZSCORESfromMATRIX($matrix,$i,$j);
#my @stores=@{$stores};

#my $scoreindex;
#if ( $scoretype eq "Z" || $scoretype eq "z") { $scoreindex = 4 } 
#elsif ( $scoretype eq "M" || $scoretype eq "m") { $scoreindex = 1 }
#else { print "ERROR: unknown --score type\n" ; exit; };

#$outfile=$matrixBed;
#$outfile=~s/.mat//;
#$outfile= $outfile . ".d" . $bflank . ".u" . $aflank . ".b" . $binwidth . ".s". $scoretype . ".zcr";
#print "$outfile\n";
#open ZWRITE,">$outfile";

#for (my $n=-$bflank;$n<=$aflank;$n+=$binwidth) {
#	my $i = (( $n - -$bflank ) / $binwidth) + 1  ;
#	print ZWRITE "$n", "\t",$$stores[$i][$scoreindex],"\n";
#}
#exit;

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
  -b, --bin        Size of each bin in bp 
  -d, --down       Window downstream of the reference coordinates. 
  -u, --up         Window upstream of the reference coordinates. 
  -s, --score      output mean(m) or z-scores(z) (default z)

Examples: 
\$ perl $fn -mat MATRIXBED -bin BIN -down DOWNSTREM -up UPSTREAM

DES
exit;
}

