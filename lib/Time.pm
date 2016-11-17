#!/usr/bin/perl -w

package Time;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw( PRINTTIME convert_time );

sub PRINTTIME {
        my $msg = shift @_;
        my $start_time = shift @_;
        my $time_tag = shift @_;
        my $time = time - $time_tag ;
        my $total_time = time - $start_time ;
        $time = convert_time($time);
        $total_time = convert_time($total_time);
        print ("** Total time elapsed: $total_time;\t $msg in $time\n\n");
        $time_tag = time;
	return($time_tag);
}

sub convert_time {
	my $time = shift;
	my $days = int($time / 86400);
	$time -= ($days * 86400);
	my $hours = int($time / 3600);
	$time -= ($hours * 3600);
	my $minutes = int($time / 60);
	my $seconds = $time % 60;

	$days = $days < 1 ? '' : $days .'d ';
	$hours = $hours < 1 ? '' : $hours .'h ';
	$minutes = $minutes < 1 ? '' : $minutes . 'm ';
	$time = $days . $hours . $minutes . $seconds . 's';
	return $time;
}
1;

