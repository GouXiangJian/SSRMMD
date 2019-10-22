#!/usr/bin/env perl

#author : Xiang-jian Gou
#date   : 2019-03-30
#email  : xjgou@stu.sicau.edu.cn

#load modules
use strict;
use warnings;
use 5.010;
use File::Basename qw/basename/;
use File::Spec;
use Storable qw/store retrieve/;
use threads;
use Getopt::Long;

#clear buffer
$| = 1;

#record version information
my $VERSION = 'SSRMMD v1.0';

#set default options
my ($fasta1_name, $fasta2_name, $genome1_name, $genome2_name, $version, $help);
my $motifs      = '1=10,2=7,3=6,4=5,5=4,6=4';
my $fs_len      = 150;
my $poly        = 0;
my $inter_file  = 0;
my $stat_file   = 0;
my $threads_num = 1;

#get options from command line
GetOptions(
        'fasta1|f1=s'    =>   \$fasta1_name,
        'fasta2|f2=s'    =>   \$fasta2_name,
        'genome1|g1=s'   =>   \$genome1_name,
        'genome2|g2=s'   =>   \$genome2_name,
        'motifs=s'       =>   \$motifs,
        'length=i'       =>   \$fs_len,
        'poly=i'         =>   \$poly,
        'detail=i'       =>   \$inter_file,
        'stat=i'         =>   \$stat_file,
        'thread=i'       =>   \$threads_num,
        'version+'       =>   \$version,
        'help+'          =>   \$help,
);

#describe program information
my $usage = <<__GUIDE__;
###############################################################################
Name:
  SSRMMD - Simple Sequence Repeat Molecular Marker Developer

Function:
  Find out SSR loci and candidate polymorphic SSR molecular marker

Usage:
  perl SSRMMD.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options: 
  -f1|-fasta1  <str> : FASTA file use to find SSR loci (must be provided)
  -f2|-fasta2  <str> : FASTA file (essential if need to find polymorphic loci)
  -g1|-genome1 <str> : genome1 file of fasta1 (default: fasta1 file)
  -g2|-genome2 <str> : genome2 file of fasta2 (default: fasta2 file)
  -m |-motifs  <str> : threshold of motifs (default: 1=10,2=7,3=6,4=5,5=4,6=4)
                       [ left  of equal : length of motif          ]
                       [ right of equal : minimum number of repeat ]
  -l |-length  <int> : length of SSR flanking sequence (default: 150)
  -p |-poly    <int> : 0 = find SSR loci, 1 = find polymorphic loci (default:0)
  -d |-detail  <int> : 1 = output intermediate file, 0 = not output (default:0)
  -s |-stat    <int> : 1 = output statistical  file, 0 = not output (default:0)
  -t |-thread  <int> : the number of threads for the program to run (default:1)
  -v |-version       : show the version information
  -h |-help          : show the help information
###############################################################################

__GUIDE__

#check the options
die "$VERSION\n" if $version;
die $usage if $help;
die "Error: option '-f1' must be provided !\n$usage" if ! $fasta1_name;
die "Error: option '-f2' must be provided !\n" if $poly == 1 and ! $fasta2_name;
die "Error: option '-l'  must be an integer bigger than 0 !\n" if $fs_len <= 0 or $fs_len =~ /\.\d*[1-9]+/;
die "Error: option '-p'  just be set to 0 or 1 !\n" if $poly != 0 and $poly != 1;
die "Error: option '-d'  just be set to 0 or 1 !\n" if $inter_file != 0 and $inter_file != 1;
die "Error: option '-s'  just be set to 0 or 1 !\n" if $stat_file != 0 and $stat_file != 1;
die "Error: option '-t'  must be an integer bigger than 0 !\n" if $threads_num <= 0 or $threads_num =~ /\.\d*[1-9]+/;

#specify genome file
$genome1_name = $fasta1_name if ! $genome1_name;
$genome2_name = $fasta2_name if $fasta2_name and ! $genome2_name;

#get base name of fasta file
my $fa_base_name1 = basename $fasta1_name;
my $fa_base_name2 = basename $fasta2_name if $fasta2_name;

#========================================#
#        start main program ...          #
#========================================#

#Module1 : find SSR loci

if($poly == 0){

	my $fa_base_name = basename $fasta1_name;

	my $molen_minum = &get_motifs_info($motifs);

	my $id_seq = &deal_fa_file($fasta1_name);

	my ($tmp_dir_name1, $final_threads_num, $id_len) = &split_id_seq_to_diff_file($id_seq, $threads_num, $fa_base_name);

	undef $id_seq;

	my $all_SSRs = &multithreads_find_SSR_and_info($molen_minum, $tmp_dir_name1, $final_threads_num, $fs_len);

	&some_simple_statistics($fa_base_name, $all_SSRs, $id_len, $molen_minum) if $stat_file;

	&generate_SSRs_with_fs_stat_file($fa_base_name, '.all_SSRs_with_fs', $all_SSRs);

}

#Module2 : develop candidate polymorphic SSR

else{

	my @retrieve_files;

	foreach my $fa_file_name ($fasta1_name, $fasta2_name){

		#get base name of fasta file
		my $fa_base_name = basename $fa_file_name;

		#step1 : find all SSRs loci
		my $molen_minum = &get_motifs_info($motifs);

		my $id_seq = &deal_fa_file($fa_file_name);

		my ($tmp_dir_name1, $final_threads_num, $id_len) = &split_id_seq_to_diff_file($id_seq, $threads_num, $fa_base_name);

		my $genome_name = $fa_file_name eq $fasta1_name ? $genome1_name : $genome2_name;

		push @retrieve_files, $genome_name eq $fa_file_name ? StoreToFile($fa_base_name, 'seq', $id_seq) : '';

		undef $id_seq;

		my $all_SSRs = &multithreads_find_SSR_and_info($molen_minum, $tmp_dir_name1, $final_threads_num, $fs_len);

		&generate_all_SSRs_stat_file($fa_base_name, $all_SSRs);

		&some_simple_statistics($fa_base_name, $all_SSRs, $id_len, $molen_minum) if $stat_file;

		&generate_SSRs_with_fs_stat_file($fa_base_name, '.all_SSRs_with_fs', $all_SSRs) if $inter_file;

		push @retrieve_files, StoreToFile($fa_base_name, 'ssr', $all_SSRs);

		undef $all_SSRs;

		undef $id_len;
	}
		#step2 : keep all flanking sequence conservative SSRs
		my $all_SSRs_add_fs1 = RetrieveFromFile($retrieve_files[1]);

		my $fs_count1 = &generate_fs_list($all_SSRs_add_fs1);

		my $fs_uni_SSRs1 = &keep_fs_uniq_SSRs($fs_count1, $all_SSRs_add_fs1);

		undef $fs_count1;

		undef $all_SSRs_add_fs1;


		my $all_SSRs_add_fs2 = RetrieveFromFile($retrieve_files[3]);

		my $fs_count2 = &generate_fs_list($all_SSRs_add_fs2);

		my $fs_uni_SSRs2 = &keep_fs_uniq_SSRs($fs_count2, $all_SSRs_add_fs2);

		undef $fs_count2;

		undef $all_SSRs_add_fs2;


		my ($fs_con_SSRs1, $fs_con_SSRs2, $con_compared) = &keep_fs_cons_SSRs($fs_uni_SSRs1, $fs_uni_SSRs2);

		undef $fs_uni_SSRs1;

		undef $fs_uni_SSRs2;

		
		&generate_cons_or_uniq_SSRs_stat_file($fa_base_name1, '.con_SSRs', $fs_con_SSRs1) if $inter_file;
		
		&generate_cons_or_uniq_SSRs_stat_file($fa_base_name2, '.con_SSRs', $fs_con_SSRs2) if $inter_file;


		my $fs_con_SSRs_file1 = StoreToFile($fa_base_name1, 'con', $fs_con_SSRs1);
		
		undef $fs_con_SSRs1;


		my $fs_con_SSRs_file2 = StoreToFile($fa_base_name2, 'con', $fs_con_SSRs2);
		
		undef $fs_con_SSRs2;


		my $con_compare_file = StoreToFile("$fa_base_name1-$fa_base_name2", 'cmp', $con_compared);
		
		undef $con_compared;


		#step3 : keep all unique SSRs
		my $id_seq1 = $genome1_name eq $fasta1_name ? RetrieveFromFile($retrieve_files[0]) : &deal_fa_file($genome1_name);

		my $tmp_dir_name2_1 = &split_id_seq_values_to_diff_file($id_seq1, $threads_num, $fa_base_name1);

		undef $id_seq1;

		my $fs_uniq_SSRs1 = &multithreads_keep_unique_SSRs($tmp_dir_name2_1, $fs_con_SSRs_file1, $fs_len, $fa_base_name1);

		&generate_cons_or_uniq_SSRs_stat_file($fa_base_name1, '.uni_SSRs', $fs_uniq_SSRs1) if $inter_file;


		my $uniq_SSRs_file1 = StoreToFile($fa_base_name1, 'uni', $fs_uniq_SSRs1);

		undef $fs_uniq_SSRs1;


		my $id_seq2 = $genome2_name eq $fasta2_name ? RetrieveFromFile($retrieve_files[2]) : &deal_fa_file($genome2_name);

		my $tmp_dir_name2_2 = &split_id_seq_values_to_diff_file($id_seq2, $threads_num, $fa_base_name2);

		undef $id_seq2;

		my $fs_uniq_SSRs2 = &multithreads_keep_unique_SSRs($tmp_dir_name2_2, $fs_con_SSRs_file2, $fs_len, $fa_base_name2);
		
		&generate_cons_or_uniq_SSRs_stat_file($fa_base_name2, '.uni_SSRs', $fs_uniq_SSRs2) if $inter_file;


		$fs_uniq_SSRs1 = RetrieveFromFile($uniq_SSRs_file1);

		$con_compared = RetrieveFromFile($con_compare_file);


		#step4 : output a table of candidate polymorphic SSRs
		&final_comparison_table($fs_uniq_SSRs1, $fs_uniq_SSRs2, $con_compared, $fa_base_name1, $fa_base_name2);
}

#========================================#
#          end main program !!!          #
#========================================#

#here are all subroutines are used by this program:

#function1 : store information to tmp file, and return a file name
sub StoreToFile {
	my ($file, $suffix, $info) = @_;
	my ($sec, $min, $hour, $day, $mon, $year, undef, undef, undef) = localtime;
	$mon  += 1;
	$year += 1900;
	my $time = "$year$mon$day-$hour$min$sec";
	my $name = "$file.$time.$suffix.tmp";
	unlink $name if -e $name;
	store $info, $name;
	undef $info;
	return $name;
}

#function2 : retrieve information from tmp file, and return the information
sub RetrieveFromFile {
	my $file = shift;	
	die "no find file '$file' : $!" unless -e $file;
	my $info = retrieve $file;
	unlink $file;
	return $info;
}

#function3 : read fasta file, and store the sequence name and corresponding sequence into a hash, and return a hash reference.
sub deal_fa_file{
    my $file_name = shift;
    open IN_FA, '<', $file_name or die "can't open $file_name:$!";
    my ($id, %id_seq) = ('', ());
    while(<IN_FA>){
    	s/[\r\n]+//;
    	#next unless $_;
    	/\A>(.*)/ ? ( $id = $1 =~ s/\s/_/gr ) : ( $_ = uc, s/[^A-Z]//g, $id_seq{$id} .= $_ );
    }
    close IN_FA;
    die "the format of fasta file may be incorrect, please check the fasta file: $file_name\n" if keys %id_seq == 0;
    return \%id_seq;
}

#function4 : split hash id_seq into smaller hashes, and put these smaller hashes in temporary files.
sub split_id_seq_to_diff_file{
	my ($id_seq, $threads_num, $fa_base_name) = @_;
	my ($sec, $min, $hour, $day, $mon, $year, undef, undef, undef) = localtime;
	$mon  += 1;
	$year += 1900;
	my $time = "$year$mon$day-$hour$min$sec";
	my $tmp_dir_name1 = $fa_base_name."_$time.SSRMMD_tmp1";
	if(-e $tmp_dir_name1){
		chdir $tmp_dir_name1;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name1;
	}
	mkdir $tmp_dir_name1, 0755 or die "can't create $tmp_dir_name1 directory:$!";
	my %id_len;
	$id_len{$_} = length $id_seq->{$_} foreach keys %$id_seq;
	my $seq_sum = scalar keys %id_len;
	my $final_threads_num = ($threads_num >= $seq_sum) ? $seq_sum : $threads_num;
	my @id_to_thread = average_allocation_algorithm($final_threads_num, %id_len);
	foreach my $i (0 .. $#id_to_thread){
		my $count = $i+1;
		my $in_thread_count = 1;
		my $each_file_size = 0;
		my %seq_to_thread = ();
		foreach my $j (1 .. $#{$id_to_thread[$i]}){
			if ($each_file_size > 10_000_000) {
				my $full_name = File::Spec->catfile($tmp_dir_name1, "$count-$in_thread_count");
				store \%seq_to_thread, $full_name;
				$in_thread_count++;
				$each_file_size = 0;
				%seq_to_thread = ();
			}
			$seq_to_thread{$id_to_thread[$i][$j]} = $id_seq->{$id_to_thread[$i][$j]};
			$each_file_size += $id_len{$id_to_thread[$i][$j]};
		}
		my $full_name = File::Spec->catfile($tmp_dir_name1, "$count-$in_thread_count");
		store \%seq_to_thread, $full_name;
	}
	undef $id_seq;
	return $tmp_dir_name1, $final_threads_num, \%id_len;
}

#function5 : allocate the sequence evenly to different threads by length, and return a two-dimensional array
sub average_allocation_algorithm{
	my ($threads_num, %id_len) = @_;
	my @id = sort { $id_len{$b} <=> $id_len{$a} } keys %id_len;
	my @id_to_thread;
	if($threads_num == 1){
		$id_to_thread[0][0] = 0; #this 0 isn't important and can be changed to any value
		push @{$id_to_thread[0]}, @id;
	}
	else{
		foreach my $i (0 .. $threads_num - 1){
			$id_to_thread[$i][0] = $id_len{$id[$i]};
			$id_to_thread[$i][1] = $id[$i];
		}
		foreach my $i ($threads_num .. $#id){
			@id_to_thread = sort { $b->[0] <=> $a->[0] } @id_to_thread if $id_to_thread[-1][0] > $id_to_thread[-2][0];
			$id_to_thread[-1][0] += $id_len{$id[$i]};
			push @{$id_to_thread[-1]}, $id[$i];
		}
	}
	return @id_to_thread;
}

#function6 : store the motif length and corresponding minimum number of repeat into a hash, and return a hash reference.
sub get_motifs_info{
	my $motifs_info = shift;
	my @digit = $motifs_info =~ /(\d+)=(\d+)/g;
	die "Error: option '-m' have a problem !\n" if @digit == 0 or @digit % 2; 
	my %molen_minum = @digit;
	return \%molen_minum;
}

#function7 : judge if motif is false, 1 = false, 0 = true
sub is_false_motif{
    my $motif = shift;
    my $motif_len = length $motif;
    my $judge = 0;
    return $judge if $motif_len == 1;
    my @composite_num;
    foreach(1 .. int($motif_len/2)){
        push @composite_num, $_ if ! ($motif_len % $_);
    }
    foreach my $len (@composite_num){
        my $tmp = $motif;
        my %motif_sub;
        while($tmp){
            my $sub = substr $tmp, 0, $len;
            $motif_sub{$sub} = 1;
            substr($tmp, 0, $len) = '';
        }
        if(keys %motif_sub == 1) {
            $judge = 1;
            last;
        }
    }
    return $judge;
}

#function8 : return a hash reference that contains all SSRs of information found.
sub find_SSR_and_info{
	my ($molen_minum, $files, $fs_len) = @_;
	my %id_info;
	foreach my $file (@$files){
		my $id_seq = retrieve $file;
		foreach my $id (sort keys %$id_seq){
			#this method is suitable for short sequence, so I cut each sequence 
			#into several short sequences of 1Mb in length
			my $cut_len = 1_000_000; #the length of short sequence is 1Mb
			my @shseq_loci = ();
			my $start_loci = 0;
			foreach (0 .. int( length($id_seq->{$id}) / $cut_len )){
				my $short_seq = substr $id_seq->{$id}, $start_loci, $cut_len + 10_000; #each short sequence is supplemented by another 10,000 bp, in order to prevent SSR from being cut
				push @shseq_loci, [ $start_loci + 1, $short_seq ];
				$start_loci += $cut_len;
			}
			delete $id_seq->{$id};
			my $supple_start = $cut_len + 1;
			foreach my $i (0 .. $#shseq_loci){
				my ($seq_start, $short_seq) = @{$shseq_loci[$i]};
				my $front_len = $seq_start - 1;
				my $seq_tmp = $short_seq;
				my $len_tmp = $front_len;
				my $before_seq = $i == 0 ? '' : substr $shseq_loci[$i-1][1], 0, $cut_len;
				my $now_seq = substr $shseq_loci[$i][1], 0, $cut_len;
				my $after_seq = $i == $#shseq_loci ? '' : substr $shseq_loci[$i+1][1], 0, $cut_len;
				my $all_seq = "$before_seq$now_seq$after_seq";
				foreach my $molen (sort {$a <=> $b} keys %$molen_minum) {
					my $remain = $molen_minum->{$molen} - 1;
					while(1){
						my ($SSR, $motif) = $short_seq =~ /(([ATCG]{$molen})\g{-1}{$remain,})/;
						last if ! $SSR;
						my $SSR_len = length $SSR;
						my $rep_num = $SSR_len / $molen;
						my $start = $-[0] + 1;
						my $end = $+[0] + $front_len;
						my $cutlen = $start+$molen*($rep_num-1);
						substr($short_seq, 0, $cutlen) = '';
						$start += $front_len;
						last if $start - $seq_start + 1 >= $supple_start;
						if (! &is_false_motif($motif)) {
							my $s_start = $i == 0 ? $start : $start - ($i-1)*$cut_len;
							my $s_end   = $s_start + $SSR_len - 1;
							my $left_fs = '';
							my $left_fs_truelen = 0;
							my $right_fs = '';
							my $right_fs_truelen = 0;
							$left_fs  = ($start > $fs_len) ? substr $all_seq, $s_start-$fs_len-1, $fs_len : substr $all_seq, 0, $s_start-1;
							$left_fs_truelen  = length $left_fs;
							$right_fs = substr $all_seq, $s_end, $fs_len;
							$right_fs_truelen = length $right_fs;
							push @{$id_info{$id}}, [$id, $motif, $molen, $rep_num, $SSR_len, $start, $end, $left_fs, $left_fs_truelen, $right_fs, $right_fs_truelen];
						}
						$front_len += $cutlen;
					}
					$short_seq = $seq_tmp;
					$front_len = $len_tmp;
				}
			}
		}
	}
	my %all_SSRs;
	foreach my $id (sort keys %id_info){
		my @SSRs = sort {$a->[5] <=> $b->[5]} @{$id_info{$id}};
		push @{$all_SSRs{$id}}, $SSRs[0];
		foreach my $i (1 .. $#SSRs){
			push @{$all_SSRs{$id}}, $SSRs[$i] if abs($SSRs[$i][6]-$SSRs[$i-1][6]) > $SSRs[$i][2]-1;
		}
	}
	undef %id_info;
	return \%all_SSRs;
	#notice : Although we have taken remedial measures, there is still 
	####### : a case that a very long SSR(30,000 bp) being cut off,
	####### : but this kind of situation is quite rare so we no need to worry about it.
}

#function9 : use multithreads and call function find_SSR_and_info to find SSRs and their information, and return a hash reference.
sub multithreads_find_SSR_and_info{
	my ($molen_minum, $tmp_dir_name1, $final_threads_num, $fs_len) = @_;
	foreach my $i (1 .. $final_threads_num) {
		my @files = glob "$tmp_dir_name1/$i-*";
		my $thr = threads->create(\&find_SSR_and_info, $molen_minum, \@files, $fs_len);
	}
	my %all_SSRs;
	while(threads->list()){
    	foreach my $thr (threads->list(threads::joinable)){
        	my $SSRs = $thr->join();
			%all_SSRs = (%all_SSRs, %$SSRs);
    	}
	}
	chdir $tmp_dir_name1;
	unlink glob '*';
	chdir '..';
	rmdir $tmp_dir_name1;
	return \%all_SSRs;
}

#function10 : output a statistical file that contains all SSRs found to the current working directory.
sub generate_all_SSRs_stat_file{
	my ($fa_base_name, $all_SSRs) = @_;
	my $out_file_name = $fa_base_name.'.all_SSRs';
	open OUT,'>',$out_file_name or die "can't generate $out_file_name:$!";
	print OUT "id\tmotif\tmotif_len\trepeat_num\tsize\tstart\tend\n";
	foreach my $id (sort keys %$all_SSRs){
		foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$all_SSRs->{$id}}){
			my $out_row = join "\t", @{$SSR}[0 .. 6];
			print OUT "$out_row\n";
		}
	}
	close OUT;
	undef $all_SSRs;
}

#function11 : output some simple statistics about all SSRs loci
sub some_simple_statistics{
	my ($base, $all_SSRs, $id_len, $molen_minum) = @_;
	my %id_len = %$id_len;
    my (%id_sum, %motif_sum, %motif_num);
	$id_sum{$_} = 0 for sort keys %id_len;
	$motif_sum{$_} = 0 for sort keys %$molen_minum;
	foreach my $id (sort keys %$all_SSRs){
		foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$all_SSRs->{$id}}){
			$id_sum{$SSR->[0]}++;
			$motif_sum{$SSR->[2]}++;
        	my $area = do{
        	    if   ($SSR->[3] <   5) {1}
        	    elsif($SSR->[3] <=  7) {2}
        	    elsif($SSR->[3] <= 10) {3}
        	    elsif($SSR->[3] <= 15) {4}
        	    elsif($SSR->[3] <= 20) {5}
        	    elsif($SSR->[3] <= 25) {6}
        	    elsif($SSR->[3] <= 30) {7}
        	    elsif($SSR->[3] <= 40) {8}
        	    else                   {9}
        	};
        	$motif_num{$SSR->[1]}{$area} += 1;
        	$motif_num{$SSR->[1]}{total} += 1;
        	$motif_num{$SSR->[1]}{rep_n} += $SSR->[3];
		}
	}
    my $out_file_name = $base.'.stat';
    open OUT_ALL, '>', $out_file_name or die "can't generate $out_file_name:$!";
	print OUT_ALL "Some simple statistics about all SSRs loci\n";
	print OUT_ALL "==========================================\n\n\n";
	print OUT_ALL "1. Number of SSRs loci of per sequence\n";
	print OUT_ALL "======================================\n\n";
	print OUT_ALL "Seq_id\tSeq_length(bp)\tSSR_number\tSSR_density(No./Mb)\n";
	my ($id_total, $len_total) = (0, 0);
	foreach my $id (sort { $id_sum{$b} <=> $id_sum{$a} or $b cmp $a } keys %id_sum) {
		$id_total += $id_sum{$id};
		$len_total += $id_len{$id};
		print OUT_ALL "$id\t$id_len{$id}\t$id_sum{$id}\t";
		printf OUT_ALL "%.2f\n", ($id_sum{$id}*1000*1000)/$id_len{$id};
	}
	print OUT_ALL "total\t$len_total\t$id_total\t";
	printf OUT_ALL "%.2f\n\n\n", ($id_total*1000*1000)/$len_total;
	print OUT_ALL "2. Number of SSRs loci of per length of motif\n";
	print OUT_ALL "===================================\n\n";
	print OUT_ALL "Motif_length(bp)\tSSR_number\tPercentage(%)\n";
	my $motif_total = 0;
	$motif_total += $motif_sum{$_} for sort keys %motif_sum;
	foreach my $motif_len (sort {$a <=> $b} keys %motif_sum) {
		print OUT_ALL "$motif_len\t$motif_sum{$motif_len}\t";
		printf OUT_ALL "%.2f\n", ($motif_sum{$motif_len}*100)/$motif_total;
	}
	print OUT_ALL "total\t$motif_total\t100.00\n\n\n";
	print OUT_ALL "3. Number of SSRs in different number of repeat in each motif\n";
	print OUT_ALL "=============================================================\n\n";
    print OUT_ALL "Motifs\t<5\t5-7\t8-10\t11-15\t16-20\t21-25\t26-30\t31-40\t>40\tTotal\tAverage repeat number\tAverage repeat length(bp)\n";
    foreach my $motif (sort { $motif_num{$b}{total} <=> $motif_num{$a}{total} } keys %motif_num){
        print  OUT_ALL $motif, "\t";
        print  OUT_ALL defined $motif_num{$motif}{$_} ? $motif_num{$motif}{$_} : 0, "\t" foreach 1 .. 9;
        print  OUT_ALL $motif_num{$motif}{total}, "\t";
        printf OUT_ALL "%.2f\t%.2f\n", $motif_num{$motif}{rep_n}/$motif_num{$motif}{total}, ($motif_num{$motif}{rep_n}/$motif_num{$motif}{total})*length($motif);
    }
    close OUT_ALL;
	undef $all_SSRs; undef %id_len; undef $id_len; undef %id_sum; undef %motif_sum; undef %motif_num;
}

#function12 : output a statistical file that contains SSRs(add flanking sequences) information to the current working directory.
sub generate_SSRs_with_fs_stat_file{
	my ($fa_base_name, $out_file_suffix, $SSRs_info) = @_;
	my $out_file_name = $fa_base_name.$out_file_suffix;
	open OUT,'>',$out_file_name or die "can't generate $out_file_name:$!";
	print OUT "id\tmotif\tmotif_len\trepeat_num\tsize\tstart\tend\tleft_fs\tleft_fs_len\tright_fs\tright_fs_len\n";
	foreach my $id (sort keys %$SSRs_info){
		foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$SSRs_info->{$id}}){
			my $out_row = join "\t", @$SSR;
			print OUT "$out_row\n";
		}
	}
	close OUT;
	undef $SSRs_info;
}

#function13 : return a hash reference that key is the flanking sequence, value is the frequency of existence.
sub generate_fs_list{
	my $all_SSRs_add_fs = shift;
	my %fs_count;
	foreach my $id (sort keys %$all_SSRs_add_fs){
        foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$all_SSRs_add_fs->{$id}}){
			my $fs = $SSR->[7].'-'.$SSR->[9]; #connect the left flanking sequence and right flanking sequence together, and separate them by using '-'
			$fs_count{$fs}++;
        }
    }
	return \%fs_count;
	undef $all_SSRs_add_fs;
}

#function14 : return a hash reference that contains the unique SSRs of flanking sequences and their information.(compared with flanking sequences of all SSRs)
sub keep_fs_uniq_SSRs{
	my ($fs_count, $all_SSRs_add_fs) = @_;
	my %fs_uni_SSRs;
	foreach my $id (sort keys %$all_SSRs_add_fs){
        foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$all_SSRs_add_fs->{$id}}){
			my $fs = $SSR->[7].'-'.$SSR->[9];
			$fs_uni_SSRs{$fs} = [ @{$SSR}[0 .. 6] ] if $fs_count->{$fs} == 1;
        }
    }
	undef $all_SSRs_add_fs;
	undef $fs_count;
	return \%fs_uni_SSRs;
}

#function15 : return three hash references that contain conservative SSRs in each file and conservative SSRs statistical table, respectively.
sub keep_fs_cons_SSRs{
	my ($fs_uni_SSRs1, $fs_uni_SSRs2) = @_;
	my (%fs_con_SSRs1, %fs_con_SSRs2, %con_compared);
	foreach my $fs1 (sort keys %$fs_uni_SSRs1){
		if(exists $fs_uni_SSRs2->{$fs1} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$fs1}[1]){ #the flanking sequence is found in every file, and SSR motif is the same
			$fs_con_SSRs1{$fs1} = $fs_uni_SSRs1->{$fs1};
			$fs_con_SSRs2{$fs1} = $fs_uni_SSRs2->{$fs1};
			$con_compared{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$fs1}} ];
		}
	}
	undef $fs_uni_SSRs1;
	undef $fs_uni_SSRs2;
	return \%fs_con_SSRs1, \%fs_con_SSRs2, \%con_compared;
}

#function16 : output a statistical file that contains conservative or unique SSRs to the current working directory.
sub generate_cons_or_uniq_SSRs_stat_file{
	my ($fa_base_name, $out_file_suffix, $SSRs_set) = @_;
	my $out_file_name = $fa_base_name.$out_file_suffix;
	open OUT,'>',$out_file_name or die "can't generate $out_file_name:$!";
	print OUT "id\tmotif\tmotif_len\trepeat_num\tsize\tstart\tend\tleft_fs\tleft_fs_len\tright_fs\tright_fs_len\n";
	foreach my $fs (sort {$SSRs_set->{$a}[0] cmp $SSRs_set->{$b}[0] or $SSRs_set->{$a}[5] <=> $SSRs_set->{$b}[5]} keys %$SSRs_set){
		my ($left_fs, $right_fs) = split /-/,$fs;
		my $left_fs_len  = length $left_fs;
		my $right_fs_len = length $right_fs;
		my $out_row = join "\t", (@{$SSRs_set->{$fs}}, $left_fs, $left_fs_len, $right_fs, $right_fs_len);
		print OUT "$out_row\n";
	}
	close OUT;
	undef $SSRs_set;
}

#function17 : the values of hash id_seq are connected together and stored in a scalar, and this scalar are divided equally into smaller scalars, finally, put these smaller scalars in temporary files, and return the directory name created.
sub split_id_seq_values_to_diff_file{
	my ($id_seq, $threads_num, $fa_base_name) = @_;
	my ($sec, $min, $hour, $day, $mon, $year, undef, undef, undef) = localtime;
	$mon  += 1;
	$year += 1900;
	my $time = "$year$mon$day-$hour$min$sec";
	my $tmp_dir_name2 = $fa_base_name."_$time.SSRMMD_tmp2";
	if(-e $tmp_dir_name2){
		chdir $tmp_dir_name2;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name2;
	}
	mkdir $tmp_dir_name2, 0755 or die "can't create $tmp_dir_name2 directory:$!";
    my $seqs = '';
	foreach my $id (sort keys %$id_seq) {
		$seqs .= $id_seq->{$id}.'---';
		delete $id_seq->{$id};
	}
	undef $id_seq;
	my $shseq_len =  int((length $seqs) / $threads_num)+1;
	my ($short_seq, $count) = ('',1);
    while($seqs){
	    $short_seq = substr $seqs, 0, $shseq_len; #notice :SSR at the boundary may be cut off, but the impact is very small
        my $new_name = $fa_base_name.'.tmp'.$count;
        my $out_file_name = File::Spec->catfile($tmp_dir_name2,$new_name);
		store \$short_seq, $out_file_name;
        $count++;
	    substr($seqs, 0, $shseq_len) = '';
    }
	undef $short_seq;
	undef $seqs;
	return $tmp_dir_name2;
}

#function18 : return a array reference that contains match count of each flanking sequence compared with genome sequence.
sub get_fs_match_count{
    my ($file, $new_fs_set, $fs_len) = @_;
	my $short_seq = retrieve $file;
	my $shseq_len = length $$short_seq;
	my @match_count = (0) x @$new_fs_set;
	foreach my $first_start_loci (0 .. ($fs_len-1)){ #the number of hash %seq_count generated is the same as the length of the flanking sequence
		my %seq_count = ();
		my $start_loci = $first_start_loci;
		foreach my $num (1 .. int( ($shseq_len - $first_start_loci) / $fs_len )){ #the number of cuts when the genome sequence is cut to the same length as the flanking sequence
			my $seq = substr $$short_seq, $start_loci, $fs_len;
			$seq_count{$seq}++;
			$start_loci += $fs_len;
		}
		my $sub = -1;
		foreach my $i (0 .. $#{$new_fs_set}){ #each flanking sequence is aligned with the genome sequence
			$sub++;
			my $count = 0;
			foreach my $each_fs (@{$new_fs_set->[$i]}){ #flanking sequence is divided into left and right
				next unless $seq_count{$each_fs};
				if($seq_count{$each_fs} > 1){ #flanking sequence on either side is matched to more than one position
					$count = 3;
					last;
				}
				$count += 1;
			}
			$match_count[$sub] += $count;
		}
	}
	undef $new_fs_set;
	undef $short_seq;
	return \@match_count;
	#notice : By doing this, the flank sequence(140bp) that is not enough to set the value(150bp) will default
	####### : to mismatch. Obviously, this kind of situation rarely happen. However, we have taken remedial 
	####### : measure: if the final match count is less than or equal to 2, the flanking sequence will be regarded as unique.
}

#function19 : use multithreads and call function get_fs_match_count to keep the unique SSRs(the flanking sequence is unique) in the genome, and return a hash reference.
sub multithreads_keep_unique_SSRs{
	my ($tmp_dir_name2, $fs_con_SSRs_file, $fs_len, $fa_base_name) = @_;
	die "no find file '$fs_con_SSRs_file' : $!" unless -e $fs_con_SSRs_file;
	my $fs_con_SSRs = retrieve $fs_con_SSRs_file;
	my $fs_set = [ sort {$fs_con_SSRs->{$a}[0] cmp $fs_con_SSRs->{$b}[0] or $fs_con_SSRs->{$a}[5] <=> $fs_con_SSRs->{$b}[5]} keys %$fs_con_SSRs ];
	undef $fs_con_SSRs;
	my $new_fs_set;
	foreach my $fs (@$fs_set){
		my @fs = split /-/, $fs;
		push @{$new_fs_set}, [@fs];
	}
	unless(defined $new_fs_set){
		chdir $tmp_dir_name2;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name2;
		unlink $fs_con_SSRs_file;
		goto LOOP;
	}
	undef $fs_set;
	foreach my $file (glob "$tmp_dir_name2/*"){
		#my $short_seq = retrieve $file;
		my $thr = threads->create(\&get_fs_match_count, $file, $new_fs_set, $fs_len);
	}
	undef $new_fs_set;
	my @all_count = ();
	while(threads->list()){
    	foreach my $thr (threads->list(threads::joinable)){
        	my $match_count = $thr->join();
			push @all_count, $match_count;
    	}
	}
	chdir $tmp_dir_name2;
	unlink glob '*';
	chdir '..';
	rmdir $tmp_dir_name2;
	my @match_count = ();
	foreach my $i (0 .. $#{$all_count[0]}){
		my $count = 0;
		$count += $all_count[$_][$i] foreach 0 .. $#all_count;
		push @match_count, $count;
	}
	undef @all_count;
    my %fs_uni_SSRs;
	my $sub = 0;
	$fs_con_SSRs = RetrieveFromFile($fs_con_SSRs_file);
	foreach my $fs (sort {$fs_con_SSRs->{$a}[0] cmp $fs_con_SSRs->{$b}[0] or $fs_con_SSRs->{$a}[5] <=> $fs_con_SSRs->{$b}[5]} keys %$fs_con_SSRs){
		$fs_uni_SSRs{$fs} = $fs_con_SSRs->{$fs} if $match_count[$sub] <= 2; #SSR with flanking sequence matching count of 0, 1, 2 is retained
		$sub++;
	}
	undef @match_count;
	undef $fs_con_SSRs;
	LOOP:
	return \%fs_uni_SSRs;
}

#function20 : output a statistical table that contains polymorphic SSRs candidate to the current working directory.
sub final_comparison_table{
	my ($fs_uni_SSRs1, $fs_uni_SSRs2, $con_compared, $fa_base_name1, $fa_base_name2) = @_;
	my $out_file_name = $fa_base_name1.'-and-'.$fa_base_name2.'.compare';
	open OUT,'>',$out_file_name or die "can't generate $out_file_name:$!";
	print OUT "file1_id\tfile1_motif\tfile1_repeat_num\tfile1_start\tfile1_end\tfile2_id\tfile2_motif\tfile2_repeat_num\tfile2_start\tfile2_end\tleft_fs\tleft_fs_len\tright_fs\tright_fs_len\tpolymorphism\n";
	foreach my $fs (sort {$con_compared->{$a}[0] cmp $con_compared->{$b}[0] or $con_compared->{$a}[5] <=> $con_compared->{$b}[5]} keys %$con_compared){
		if(exists $fs_uni_SSRs1->{$fs} and exists $fs_uni_SSRs2->{$fs}){
			my ($left_fs, $right_fs) = split /-/, $fs;
			my $left_fs_len  = length $left_fs;
			my $right_fs_len = length $right_fs;
			my $judge = $con_compared->{$fs}[3] != $con_compared->{$fs}[10] ? 'yes' : 'no';
			my $out_row = join "\t", (@{$con_compared->{$fs}}[0,1,3,5,6, 7,8,10,12,13], $left_fs, $left_fs_len, $right_fs, $right_fs_len, $judge);
			print OUT "$out_row\n";
		}
	}
	close OUT;
	undef $fs_uni_SSRs1; undef $fs_uni_SSRs2; undef $con_compared;
}
