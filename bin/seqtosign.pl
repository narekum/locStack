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

#use lib "/home/naren/build/punarvasu/lib" ;
use lib "$dir/../lib" ;
#use strict ;
use LWP::Simple;
use Getopt::Long;
use Parallel::ForkManager;
use GenomeInfo;
use Chatta;
use Time;

use SeqToProfile;

print '                                                                  

            ____                                                         
           / __ \ __  __ ____   ____ _ _____ _   __ ____ _ _____ __  __  
          / /_/ // / / // __ \ / __ `// ___/| | / // __ `// ___// / / /  
         / ____// /_/ // / / // /_/ // /    | |/ // /_/ /(__  )/ /_/ /   
        /_/     \__,_//_/ /_/ \__,_//_/     |___/ \__,_//____/ \__,_/    
        **            Meta-analysis of genomic features           **     


';

my $cmd=$0." ".join(" ",@ARGV); ### command line copy

my $help;
my $exp;
my $summitfile;
#my $frequencyfile;
my $outfile;
my $bflank=100;                                 
my $aflank=100;                                
my $nproc=2;                                 
my $genome="hg19";                       

$time_tag=$start_time=time;

GetOptions ('h|help'=>\$help,                     # --help         : print this help
            "s|seq=s" => \$exp,                   # --aln          : bed file of aligned reads  (required)
            "r|ref=s" => \$summitfile,            # --ref          : bed file on reference reads (required)
#            "f|freqfile=s" => \$frequencyfile,    # --freqfile     : cut frequencies of nucleotides
            "o|outfile=s" => \$outfile,           # --outfile      : outfile file name (required)
            "g|genome=s" => \$genome,             # --genome       : genome assembly name eq "hg19", mm9 etc.
            "d|down=i" => \$bflank,               # --down         : Window downstream of the reference coordinates (default 1000)
            "u|up=i" => \$aflank,                 # --up           : Window upstream of the reference coordinates (default 1000)
#            "t|thread=i" => \$nproc,              # --thread       : number of processors to be used (default 2)
           )
or &PRINTHELP($fn); 


if(defined $help || !$exp || !$genome || !$outfile || !$summitfile ) {
	&PRINTHELP($fn);
	exit;
}

####################################################################

my @options=( "$cmd",
           "--seq            $exp",
           "--ref            $summitfile",
           "--outfile        $outfile",
           "--genome         $genome",
           "--down           $bflank",
           "--up             $aflank",
#           "--freqfile       $frequencyfile",
);

print "** Running $fn with the following options\n\n" ;
print "      $_\n" foreach @options ;
print "\n";


$time_tag=Time::PRINTTIME("Options loaded in ",$start_time,$time_tag);

($summitsFwd,$summitsRev,$totalfwdsummits,$totalrevsummits)=Chatta::READSUMMITS($summitfile);

$time_tag=Time::PRINTTIME("Loaded reference coordinates ",$start_time,$time_tag);

@bincontainer = Chatta::BINS(1,$bflank,$aflank);

#%freqhash=SeqToProfile::READFREQUENCIES($frequencyfile);
#print "$_\t$freqhash{$_} \n" foreach keys %freqhash ;


#$seq="ATGCGTCG";
#$prefbed=SeqToProfile::SeqToPref($seq,\%freqhash,6);
#print "$_\t" foreach @$prefbed ; print "\n";

#my $matrixprefbed=&SeqToProfile::BedToPref($summitsFwd,$bflank,$aflank,6,$genome,$exp,\%freqhash);
my ($matrixprefbed_plus_same,$matrixprefbed_plus_oppo)=&SeqToProfile::BedToPref($summitsFwd,$bflank,$aflank,6,$genome,$exp);
my ($matrixprefbed_minus_same,$matrixprefbed_minus_oppo)=&SeqToProfile::BedToPref($summitsRev,$bflank,$aflank,6,$genome,$exp);

open SAME,">$outfile.hexa.same.mat";
print SAME "$_\n" foreach @$matrixprefbed_plus_same,@$matrixprefbed_minus_same;

open OPPO,">$outfile.hexa.oppo.mat";
print OPPO "$_\n" foreach @$matrixprefbed_plus_oppo, @$matrixprefbed_minus_oppo ;

open BOTH,">$outfile.hexa.both.mat";
print BOTH "$_\n" foreach @$matrixprefbed_plus_same,@$matrixprefbed_minus_same,@$matrixprefbed_plus_oppo, @$matrixprefbed_minus_oppo ;

#print "$_\n" foreach @$matrixprefbed ;


#print "      Begining Chromosome wise distribution of tag counts in $totalfwdsummits coordinates on plus strand in ref file $summitfile\n" ;
#print "      From $bflank to $aflank with binwidth of $binwidth; ", scalar @bincontainer , " columns.\n\n" ;

#my $FwdDist =Chatta::RETURNBINCOUNTSALLCHR ($summitsFwd,$alignedreadhash,$control,\@bincontainer,$binwidth,$extend_length,$normalize,$totalExpReads,$totalExpReads,$nproc,$genome);

#$time_tag=Time::PRINTTIME("Distribution of $totalfwdsummits coordinates on plus strand in ref file $summitfile completed ",$start_time,$time_tag);

#print "      Begining Chromosome wise distribution of tag counts in $totalrevsummits coordinates on minus strand in ref file $summitfile\n" ;
#print "      From $bflank to $aflank with binwidth of $binwidth; ", scalar @bincontainer , " columns.\n\n" ;

#my $RevDist =Chatta::RETURNBINCOUNTSALLCHR ($summitsRev,$alignedreadhash,$control,\@bincontainer,$binwidth,$extend_length,$normalize,$totalExpReads,$totalExpReads,$nproc,$genome);

#$time_tag=Time::PRINTTIME("Distribution of $totalrevsummits coordinates on minus strand in ref file $summitfile completed ",$start_time,$time_tag);

#open BOTH,">$outfile.both.mat";
#open PLUS,">$outfile.same.mat";
#open MINUS,">$outfile.oppo.mat";

#my ($FwdBothDist,$FwdPlusDist,$FwdMinusDist)=@$FwdDist;
#my ($RevBothDist,$RevPlusDist,$RevMinusDist)=@$RevDist;

#print PLUS "$_\n" foreach ( @$FwdPlusDist , @$RevMinusDist );
#print MINUS "$_\n" foreach ( @$FwdMinusDist , @$RevPlusDist ) ;
#print BOTH "$_\n" foreach ( @$FwdBothDist , @$RevBothDist );

#$time_tag=Time::PRINTTIME("Distributions printed ",$start_time,$time_tag);

#print "Total time taken in processing ", Time::convert_time($time_tag- $start_time), "; Program ending...OK.\n\n";



#system "./zscores.pl $outfile.$ori.both.mat > $outfile.$ori.both.mat.dat" ;
#system "./zscores.pl $outfile.$ori.plus.mat > $outfile.$ori.plus.mat.dat" ;
#system "./zscores.pl $outfile.$ori.minus.mat > $outfile.$ori.minus.mat.dat" ;





sub PRINTHELP {
	my ($fn)=@_;

        print <<DES;
Usage: $fn < -a ALIGNEDREADS -r REFERENCECOORDS -o NAME -g GENOME > [Options]

Options:
  -h, --help       Show this help message.
  -s, --seq        fasta file of genomic sequences  (REQUIRED)
  -r, --ref        Bed file on reference coordinates (REQUIRED)
  -o, --outfile    outputfile (REQUIRED)
  -g, --genome     genome assembly name eq "hg19", mm9 etc.
                   Available genomes:dm2,dm3,hg18,hg19,mm8,mm9,pombe,rn4 and sacCer1
  -d, --down       Window downstream of the reference coordinates. 
                   (default 1000)
  -u, --up         Window upstream of the reference coordinates. 
                   (default 1000)

Examples: 
\$ perl $fn --seq genomicseqfile --ref reference.bed --outfile outfileprefix  --genome hg19 --down 100 --up 100 
\\or\\
\$ perl $fn  -s genomicseqfile -r reference.bed -o outfileprefix -g hg19 -d 100 -u 100

DES
exit;
}
