#!/usr/bin/env perl

#usage: perl SSR_stater.pl <SSR_file>
#note : the SSR_stater.pl only accept the file with suffix .all_SSRs outputed from SSR_finder.pl

use strict;
use warnings;
use File::Basename qw/basename/;

if($ARGV[0] =~ /.all_SSRs$/){
    &stat_SSR_num_in_each_motif($ARGV[0]);
}
else{
    die "The input file may not meet the requirements\n";
}

#function : output the number of SSRs in different repeat classes in each motif
sub stat_SSR_num_in_each_motif{
    my $in_file_name = shift;
    open IN_ALL,'<',$in_file_name or die "can't open $in_file_name:$!";
    my %motif_num;
    $_ = <IN_ALL>;
    while(<IN_ALL>){
        my @row = split;
        my $area = do{
            if($row[3] <  5)     {1}
            elsif($row[3] <=  7) {2}
            elsif($row[3] <= 10) {3}
            elsif($row[3] <= 15) {4}
            elsif($row[3] <= 20) {5}
            elsif($row[3] <= 25) {6}
            elsif($row[3] <= 30) {7}
            elsif($row[3] <= 40) {8}
            else                 {9}
        };
        $motif_num{$row[1]}{$area} += 1;
        $motif_num{$row[1]}{total} += 1;
        $motif_num{$row[1]}{rep_n} += $row[3];
    }
    close IN_ALL;
    my $base_name = basename $in_file_name;
    my $out_file_name = $base_name.'.stat';
    open OUT_ALL, '>', $out_file_name or die "can't generate $out_file_name:$!";
    print OUT_ALL "Motifs\t<5\t5~7\t8~10\t11~15\t16~20\t21~25\t26~30\t31~40\t>40\tTotal\tAverage repeat number\tAverage repeat length(bp)\n";
    foreach my $motif (sort { $motif_num{$b}{total} <=> $motif_num{$a}{total} } keys %motif_num){
        print  OUT_ALL $motif, "\t";
        print  OUT_ALL defined $motif_num{$motif}{$_} ? $motif_num{$motif}{$_} : 0, "\t" foreach 1 .. 9;
        print  OUT_ALL $motif_num{$motif}{total}, "\t";
        printf OUT_ALL "%.2f\t%.2f\n", $motif_num{$motif}{rep_n}/$motif_num{$motif}{total}, ($motif_num{$motif}{rep_n}/$motif_num{$motif}{total})*length($motif);
    }
    close OUT_ALL;
}
