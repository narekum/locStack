#!/usr/bin/perl -w

package CreateIndex;
use Fcntl;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw( BUILD_INDEX LINE_WITH_INDEX CHR_OFFSET CREATE_OFFSET_INDEX);

sub RETURN_LINE {
	my ($filename,$line_number)=@_;

	open(ORIG, "< $filename") or die "Can't open $filename for reading: $!";

	my $indexname = "$filename.index";
	my $indexoffset = "$filename.offset";
	sysopen(IDX, $indexname, O_CREAT|O_RDWR,0666) or die "Can't open $indexname for read/write: $!";	

	if ( -z $indexname ) {
		&CREATE_OFFSET_INDEX($filename);
	}

	my $line = LINE_WITH_INDEX(*ORIG, *IDX, $line_number);
	return ($line);
}

sub CREATE_OFFSET_INDEX {
	my ($filename)=@_;

	open(ORIG, "< $filename")
        or die "Can't open $filename for reading: $!";

	my $indexname = "$filename.index";
	sysopen(IDX, $indexname, O_CREAT|O_RDWR,0666) or die "Can't open $indexname for read/write: $!";

	&BUILD_INDEX(*ORIG, *IDX) if -z $indexname;  # XXX: race unless lock
	&CHR_OFFSET($filename);
}


sub CHR_OFFSET {
        my ($bedfile)=@_;

        my $offset_index="$bedfile.offset";
        if ( -e $offset_index) {
                return()
        } else {
                print "  $offset_index not present -- creating it now..\n";
        }

        my $i=0;
        my $j=0;
        my $old_chr="";
        my $chr;
        my @chr_offset=();
        open OPEN,$bedfile or die "Unable to open $bedfile:$!\n";
        while (<OPEN>){
                ($chr)=split(/\t/,$_,0);
                if ($chr ne $old_chr) {
                        push @chr_offset, "$old_chr\t$j\t$i" if $old_chr ne "";
                        $j=$i+1;
                }
                $i++;
                #print "$chr\n";
                $old_chr=$chr;
        }
        push @chr_offset, "$old_chr\t$j\t$i";

        open WRITE,">$offset_index";
        print WRITE "$_\n" foreach @chr_offset ;
        return ;

}

sub BUILD_INDEX {
    my $data_file  = shift;
    my $index_file = shift;
    my $offset     = 0;
	print "  Building Index.....";
    while (<$data_file>) {
        print $index_file pack("N", $offset);
        $offset = tell($data_file);
    }
}

sub LINE_WITH_INDEX {
    my $data_file   = shift;
    my $index_file  = shift;
    my $line_number = shift;

    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file

    $size = length(pack("N", 0));
    $i_offset = $size * ($line_number-1);
    seek($index_file, $i_offset, 0) or return;
    read($index_file, $entry, $size);
    $d_offset = unpack("N", $entry);
    seek($data_file, $d_offset, 0);
    return scalar(<$data_file>);
}
