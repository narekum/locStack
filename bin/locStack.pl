#! /usr/bin/perl 

#############################################
##  Narendra Kumar, PhD                     #
##  Epigenetics Unit                        #
##  Institute of Cancer Sciences            #
##  University of Glasgow, UK               #
##  narekum@gmail.com                       #
#############################################

BEGIN {
        use Cwd qw(realpath cwd);
        use File::Basename;
        our ($fn, $dir) = fileparse(realpath($0));
}

use lib "$dir/../lib" ;
#use strict ;
use LWP::Simple;
use Getopt::Long;
use Parallel::ForkManager;
use GenomeInfo;
use Chatta;
use Time;

print '                                                                  
           _             _____ _             _    
          | |           / ____| |           | |   
          | | ___   ___| (___ | |_ __ _  ___| | __
          | |/ _ \ / __|\___ \| __/ _` |/ __| |/ /
          | | (_) | (__ ____) | || (_| | (__|   < 
          |_|\___/ \___|_____/ \__\__,_|\___|_|\_\
             Meta analysis of genomic features                                        

';

my $cmd=$0." ".join(" ",@ARGV); ### command line copy

my $binwidth=100;                                 
my $bflank=1000;                                 
my $aflank=1000;                                
my $extend_length=150;                         
my $control="NULL";                           
my $nproc=2;                                 
my $rscore="n";                             
my $dcol=5;                                
my $normalize=0;                          
my $genome="hg19";                       

$time_tag=$start_time=time;

GetOptions ('h|help'=>\$help,                     # --help         : print this help
            "a|aln=s" => \$exp,                   # --aln          : bed file of aligned reads  (required)
            "r|ref=s" => \$summitfile,            # --ref          : bed file on reference reads (required)
            "o|outfile=s" => \$outfile,           # --outfile      : outfile file name (required)
            "g|genome=s" => \$genome,             # --genome       : genome assembly name eq "hg19", mm9 etc.
            "e|ext=i" => \$extend_length,         # --ext          : how much the reads should be extended (default 100)
            "b|bin=i" => \$binwidth,              # --bin          : size of each bin in bp (default 100) 
            "d|down=i" => \$bflank,               # --down         : Window downstream of the reference coordinates (default 1000)
            "u|up=i" => \$aflank,                 # --up           : Window upstream of the reference coordinates (default 1000)
            "c|control=s" => \$control,           # --control      : Bed file of control reads (optional)
            "t|thread=i" => \$nproc,              # --thread       : number of processors to be used (default 2)
            "s|score=s" => \$rscore,              # --score        : if the read scores are to be used y=yes; n=no (default n)
            "f|field=i" => \$dcol,                # --field        : Data column in the bed file to be used if --score=y ; (default 5)
            "n|normalize=i" => \$normalize,       # --normalize    : 1000000 for reads per million (default no normalization)
           )
or &PRINTHELP($fn); 


if(defined $help || !$exp || !$genome || !$outfile || !$summitfile ) {
	&PRINTHELP($fn);
	exit;
}

####################################################################

my @options=( "$cmd",
           "--aln            $exp",
           "--ref            $summitfile",
           "--outfile        $outfile",
           "--genome         $genome",
           "--ext            $extend_length",
           "--bin            $binwidth",
           "--down           $bflank",
           "--up             $aflank",
           "--control        $control",
           "--thread         $nproc",
           "--score          $rscore",
           "--field          $dcol",
           "--normalize      $normalize",
);

print "** Running $fn with the following options\n\n" ;
print "      $_\n" foreach @options ;
print "\n";


$time_tag=Time::PRINTTIME("Options loaded in ",$start_time,$time_tag);

($alignedreadhash,$totalExpReads)=Chatta::READALIGNEDTAGS($exp,$extend_length,$rscore,$dcol);

print "     $totalExpReads tags processed from $exp file\n\n";

$time_tag=Time::PRINTTIME("Loaded aligned reads ",$start_time,$time_tag);

($summitsFwd,$summitsRev,$totalfwdsummits,$totalrevsummits)=Chatta::READSUMMITS($summitfile);

$time_tag=Time::PRINTTIME("Loaded reference coordinates ",$start_time,$time_tag);

@bincontainer = Chatta::BINS($binwidth,$bflank,$aflank);

print "      Begining Chromosome wise distribution of tag counts in $totalfwdsummits coordinates on plus strand in ref file $summitfile\n" ;
print "      From $bflank to $aflank with binwidth of $binwidth; ", scalar @bincontainer , " columns.\n\n" ;

my $FwdDist =Chatta::RETURNBINCOUNTSALLCHR ($summitsFwd,$alignedreadhash,$control,\@bincontainer,$binwidth,$extend_length,$normalize,$totalExpReads,$totalExpReads,$nproc,$genome);

$time_tag=Time::PRINTTIME("Distribution of $totalfwdsummits coordinates on plus strand in ref file $summitfile completed ",$start_time,$time_tag);

print "      Begining Chromosome wise distribution of tag counts in $totalrevsummits coordinates on minus strand in ref file $summitfile\n" ;
print "      From $bflank to $aflank with binwidth of $binwidth; ", scalar @bincontainer , " columns.\n\n" ;

my $RevDist =Chatta::RETURNBINCOUNTSALLCHR ($summitsRev,$alignedreadhash,$control,\@bincontainer,$binwidth,$extend_length,$normalize,$totalExpReads,$totalExpReads,$nproc,$genome);

$time_tag=Time::PRINTTIME("Distribution of $totalrevsummits coordinates on minus strand in ref file $summitfile completed ",$start_time,$time_tag);

open BOTH,">$outfile.both.mat";
open PLUS,">$outfile.same.mat";
open MINUS,">$outfile.oppo.mat";

my ($FwdBothDist,$FwdPlusDist,$FwdMinusDist)=@$FwdDist;
my ($RevBothDist,$RevPlusDist,$RevMinusDist)=@$RevDist;

print PLUS "$_\n" foreach ( @$FwdPlusDist , @$RevMinusDist );
print MINUS "$_\n" foreach ( @$FwdMinusDist , @$RevPlusDist ) ;
print BOTH "$_\n" foreach ( @$FwdBothDist , @$RevBothDist );

$time_tag=Time::PRINTTIME("Distributions printed ",$start_time,$time_tag);

print "Total time taken in processing ", Time::convert_time($time_tag- $start_time), "; Program ending...OK.\n\n";



#system "./zscores.pl $outfile.$ori.both.mat > $outfile.$ori.both.mat.dat" ;
#system "./zscores.pl $outfile.$ori.plus.mat > $outfile.$ori.plus.mat.dat" ;
#system "./zscores.pl $outfile.$ori.minus.mat > $outfile.$ori.minus.mat.dat" ;





sub PRINTHELP {
	my ($fn)=@_;

        print <<DES;
Usage: $fn < -a ALIGNEDREADS -r REFERENCECOORDS -o NAME -g GENOME > [Options]

Options:
  -h, --help       Show this help message.
  -a, --aln        Bed file of aligned reads  (REQUIRED)
  -r, --ref        Bed file on reference reads (REQUIRED)
  -o, --outfile    Outout file name (REQUIRED)
  -g, --genome     genome assembly name eq "hg19", mm9 etc.
                   Available genomes:dm2,dm3,hg18,hg19,mm8,mm9,pombe,rn4 and sacCer1

  -e, --ext        How much the reads should be extended 
                   (default 100)
  -b, --bin        Size of each bin in bp 
                   (default 100)
  -d, --down       Window downstream of the reference coordinates. 
                   (default 1000)
  -u, --up         Window upstream of the reference coordinates. 
                   (default 1000)
  -c, --control    Bed file of control reads (optional)
  -t, --thread     Number of processors to be used (default 2)
  -s, --score      If the read scores are to be used y=yes; n=no (default n)
  -f, --field      Data column in the bed file to be used if --score=y ; (default 5)
  -n, --normalize  1000000 for reads per million (default no normalization)

Examples: 
\$ perl $fn --aln alignedreads.bed --ref reference.bed --o profile  -genome hg19 
\\or\\
\$ perl $fn  -a alignedreads.bed -r reference.bed -o profile -g hg19 -b 50 -d 1000 -u 1000 -e 100

DES
exit;
}

