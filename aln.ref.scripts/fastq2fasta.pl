#!/usr/bin/perl

use strict;
use warnings ;

my $sequence_file = shift or die "Usage: $0 <sequence file>\n" ;
open(FILE, "<$sequence_file") or die "Failed to open file '$sequence_file'\n$!\n";
  

while (my $id_line = <FILE>) {
    if ($id_line =~ m/^@(\S+)/) {
	my $id = $1;
	my $seq_line =<FILE>;
	chomp $seq_line;
	if ($seq_line =~ m/^([ACGTN]+)$/) {
	    my $seq = $1;
	    my $second_id_line = <FILE>;
	    if ($second_id_line =~ m/^\+/) {
		my $quality_line = <FILE>;
		if ($quality_line =~ m/^.+/) {
		    my $quality = $1;
		    
		    print ">$id\n$seq\n";
		    
		} else {
		    die "Failed to parse quality line: $quality_line\n";
		}
	    } else {
		die "Failed to parse second id line: $second_id_line\n";
	    }
	} else {
	    die "Failed to parse seq line: $seq_line\n";
	}
    } else {
	die "Failed to parse id line: $id_line\n";
    }
}  



