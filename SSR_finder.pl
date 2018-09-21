#!/usr/bin/env perl

#author : Gou Xiangjian, Shi Haoran
#date   : 2018/9/11

#load modules
use strict;
use warnings;

use 5.010;
use File::Basename qw/basename/;
use File::Spec;
use Storable qw/store retrieve/;
use threads;
use Getopt::Long;

#declare options
my ($fasta1_name, $fasta2_name, $genome1_name, $genome2_name, $search_file,
    $threads_num, $fs_len, $step, $help, $version);

#get the options
GetOptions( 'fasta1=s'    =>  \$fasta1_name,
            'fasta2=s'    =>  \$fasta2_name,
            'genome1=s'   =>  \$genome1_name,
            'genome2=s'   =>  \$genome2_name,
            'search=s'    =>  \$search_file,
            'thread=i'    =>  \$threads_num,
            'length=i'    =>  \$fs_len,
            'step=i'      =>  \$step,
            'help+'       =>  \$help,
            'version+'    =>  \$version,
);

#record version information
my $VERSION = 'SSR_finder v1.0'; 

#describe program information
my $usage = <<__GUIDE__;

####################################################################################################
Function:
  find SSR loci and polymorphic SSR candidate.

Usage:
  SSR_finder.pl -fasta1 <fasta_file1> -fasta2 <fasta_file2> -se <search_condition_file> -st <step_num> -t <threads_num>

Required Options: 
  -fasta1      <file>   : A file of fasta format used to find SSR loci.
  -se|-search  <file>   : A file used to set SSR search conditions, and its content should be like this:
                          motif_len     min_repeat_num
                          1             10
                          2             7
                          3             6
                          4             5
                          5             4
                          6             4

Optional Options:
  -fasta2      <file>   : A file of fasta format used to find SSR loci. If not, the program will only look for
                          SSR loci for fasta file1, similarly, not find polymorphic SSR candidate.
  -genome1     <file>   : A file of fasta format used to specify the reference genome of fasta1. If not, it is
                          assumed that the reference genome of fasta1 is itself.
  -genome2     <file>   : A file of fasta format used to specify the reference genome of fasta2. If not, it is
                          assumed that the reference genome of fasta2 is itself.
  -t|-thread   <int>    : An integer that sets the number of threads for the program to run. (default: 1)
  -l|-length   <int>    : An integer that sets the length of SSR flanking sequence. (default: 150)
  -st|-step    <int>    : It only be set to 1 or 2. And 1 for only find SSR loci, 2 for find polymorphic SSR
                          candidate. (default: 1)
  -h|-help              : show the help information
  -v|-version           : show the version information
####################################################################################################

__GUIDE__

#some die operation
die "$VERSION\n" if $version; #show the version information

die $usage if $help; #show the help information

die $usage unless $fasta1_name and $search_file; #lack of required options

die "the value of option '-thread' must be an integer bigger than 0\n" if (defined $threads_num) and ($threads_num <= 0 or $threads_num =~ /\.\d*[1-9]+/);

die "the value of option '-length' must be an integer bigger than 0\n" if (defined $fs_len) and ($fs_len <= 0 or $fs_len =~ /\.\d*[1-9]+/);

die "the value of option '-step' only be set to 1 or 2\n" if (defined $step) and ($step != 1 and $step != 2);

die "you may need to specify another option '-fasta2'\n" if defined $step and $step == 2 and !(defined $fasta2_name);

#set the default value
$threads_num = 1    unless defined $threads_num;
$fs_len      = 150  unless defined $fs_len;
$step        = 1    unless defined $step;

#specify the reference genome file
$genome1_name = $fasta1_name unless defined $genome1_name;
$genome2_name = $fasta2_name if !(defined $genome2_name) and defined $fasta2_name;

#get base name of fasta file
my $fa_base_name1 = basename $fasta1_name;
my $fa_base_name2 = basename $fasta2_name if defined $fasta2_name;


#start main program

if($step == 1){
		#get base name of fasta file
		my $fa_base_name = basename $fasta1_name;

		#find all SSRs loci
		my $id_seq = &deal_fa_file($fasta1_name);

		my $tmp_dir_name1 = &split_id_seq_to_diff_file($id_seq, $threads_num, $fa_base_name);

		my $molen_minum = &deal_search_file($search_file);

		my $all_SSRs = &multithreads_find_SSR_and_info($molen_minum, $tmp_dir_name1);

		&generate_all_SSRs_stat_file($fa_base_name, $all_SSRs);
}
else{
	foreach my $fa_file_name ($fasta1_name, $fasta2_name){
		#get base name of fasta file
		my $fa_base_name = basename $fa_file_name;

		#step1 : find all SSRs loci
		my $id_seq = &deal_fa_file($fa_file_name);

		my $tmp_dir_name1 = &split_id_seq_to_diff_file($id_seq, $threads_num, $fa_base_name);

		my $molen_minum = &deal_search_file($search_file);

		my $all_SSRs = &multithreads_find_SSR_and_info($molen_minum, $tmp_dir_name1);

		&generate_all_SSRs_stat_file($fa_base_name, $all_SSRs);


		#step2 : add flanking sequences to all SSRs found
		my $tmp_dir_name2 = &split_SSR_and_id_seq_to_diff_file($all_SSRs, $id_seq, $threads_num, $fa_base_name);

		my $all_SSRs_add_fs = &multithreads_add_flanking_seq($fs_len, $tmp_dir_name2);

		#if you want statistical file of all SSRs with flanking sequence, please remove the # character in front of the next line
		#&generate_SSRs_with_fs_stat_file($fa_base_name, '.all_SSRs_with_fs', $all_SSRs_add_fs);


		my $genome_name = $fa_file_name eq $fasta1_name ? $genome1_name : $genome2_name;

		if($genome_name eq $fa_file_name){

			#store sequences information for each file
			unlink $fa_base_name.'.seq-tmp' if -e $fa_base_name.'.seq-tmp';

			store $id_seq, $fa_base_name.'.seq-tmp';
		}

		#store SSRs information for each file
		unlink $fa_base_name.'.ssr-tmp' if -e $fa_base_name.'.ssr-tmp';

		store $all_SSRs_add_fs, $fa_base_name.'.ssr-tmp';
	}

		#step3 : keep all flanking sequence conservative SSRs
		die "no find file ".$fa_base_name1.".ssr-tmp:$!" unless -e $fa_base_name1.'.ssr-tmp';

		my $all_SSRs_add_fs1 = retrieve $fa_base_name1.'.ssr-tmp';

		unlink $fa_base_name1.'.ssr-tmp';

		my $fs_count1 = &generate_fs_list($all_SSRs_add_fs1);

		my $fs_uni_SSRs1 = &keep_fs_uniq_SSRs($fs_count1, $all_SSRs_add_fs1);

		undef $all_SSRs_add_fs1; #memory release


		die "no find file ".$fa_base_name2.".ssr-tmp:$!" unless -e $fa_base_name2.'.ssr-tmp';

		my $all_SSRs_add_fs2 = retrieve $fa_base_name2.'.ssr-tmp';

		unlink $fa_base_name2.'.ssr-tmp';

		my $fs_count2 = &generate_fs_list($all_SSRs_add_fs2);

		my $fs_uni_SSRs2 = &keep_fs_uniq_SSRs($fs_count2, $all_SSRs_add_fs2);

		undef $all_SSRs_add_fs2; #memory release


		my ($fs_con_SSRs1, $fs_con_SSRs2, $con_compared) = &keep_fs_cons_SSRs($fs_uni_SSRs1, $fs_uni_SSRs2);

		#if you want statistical file1 of conservative SSRs, please remove the # character in front of the next line
		#&generate_cons_or_uniq_SSRs_stat_file($fa_base_name1, '.con_SSRs', $fs_con_SSRs1);
		
		#if you want statistical file2 of conservative SSRs, please remove the # character in front of the next line
		#&generate_cons_or_uniq_SSRs_stat_file($fa_base_name2, '.con_SSRs', $fs_con_SSRs2);


		#step4 : keep all unique SSRs
		my $id_seq1;

		if($genome1_name eq $fasta1_name){

			die "no find file ".$fa_base_name1.".seq-tmp:$!" unless -e $fa_base_name1.'.seq-tmp';

			$id_seq1 = retrieve $fa_base_name1.'.seq-tmp';

			unlink $fa_base_name1.'.seq-tmp';

		}
		else{

			$id_seq1 = &deal_fa_file($genome1_name);

		}

		my $tmp_dir_name3_1 = &split_id_seq_values_to_diff_file($id_seq1, $threads_num, $fa_base_name1);

		undef $id_seq1; #memory release of fasta file1

		my $fs_uniq_SSRs1 = &multithreads_keep_unique_SSRs($tmp_dir_name3_1, $fs_con_SSRs1, $fs_len);

		#if you want statistical file1 of unique SSRs, please remove the # character in front of the next line
		#&generate_cons_or_uniq_SSRs_stat_file($fa_base_name1, '.uni_SSRs', $fs_uniq_SSRs1);


		my $id_seq2;

		if($genome2_name eq $fasta2_name){

			die "no find file ".$fa_base_name2.".seq-tmp:$!" unless -e $fa_base_name2.'.seq-tmp';

			$id_seq2 = retrieve $fa_base_name2.'.seq-tmp';

			unlink $fa_base_name2.'.seq-tmp';

		}
		else{

			$id_seq2 = &deal_fa_file($genome2_name);

		}

		my $tmp_dir_name3_2 = &split_id_seq_values_to_diff_file($id_seq2, $threads_num, $fa_base_name2);

		undef $id_seq2; #memory release of fasta file2

		my $fs_uniq_SSRs2 = &multithreads_keep_unique_SSRs($tmp_dir_name3_2, $fs_con_SSRs2, $fs_len);
		
		#if you want statistical file2 of unique SSRs, please remove the # character in front of the next line
		#&generate_cons_or_uniq_SSRs_stat_file($fa_base_name2, '.uni_SSRs', $fs_uniq_SSRs2);


		#step5 : output a table of polymorphic SSRs candidate
		&final_comparison_table($fs_uniq_SSRs1, $fs_uniq_SSRs2, $con_compared, $fa_base_name1, $fa_base_name2);
}


#some subroutine

#----------------------------------------------------------------------------------------------------#

#function : read fasta file, and store the sequence name and corresponding sequence into a hash, and return a hash reference.
sub deal_fa_file{
	my $file_name = shift;
    open IN_FA,'<',$file_name or die "can't open $file_name:$!";
	my ($id, %id_seq);
	while(<IN_FA>){
    	chomp;
    	next unless $_;
    	/\A>(.*)/ ? ($id = $1 =~ s/\s/_/gr) : ($id_seq{$id} .= uc);
	}
	close IN_FA;
	die "the format of the fasta file may be incorrect, please check the fasta file: $file_name\n" if keys %id_seq == 0;
	return \%id_seq;
}

#function : split hash id_seq into smaller hashes, and put these smaller hashes in temporary files, and return the directory name created.
sub split_id_seq_to_diff_file{
	my ($id_seq, $threads_num, $fa_base_name) = @_;
	my $tmp_dir_name1 = $fa_base_name.'.find_SSR_tmp1';
	if(-e $tmp_dir_name1){
		chdir $tmp_dir_name1;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name1;
	}
	mkdir $tmp_dir_name1, 0755 or die "can't create $tmp_dir_name1 directory:$!";
	my %id_len;
	$id_len{$_} = length $id_seq->{$_}  foreach keys %$id_seq;
	my @id = sort { $id_len{$b} <=> $id_len{$a} } keys %id_len;
	my $file_num = ($threads_num >= @id) ? @id : $threads_num;
	my @id_to_thread;
	if($file_num == 1){
		$id_to_thread[0][0] = 0; #this 0 isn't important and can be changed to any value
		push @{$id_to_thread[0]}, @id;
	}
	else{
		foreach my $i (0 .. $file_num - 1){
			$id_to_thread[$i][0] = $id_len{$id[$i]};
			$id_to_thread[$i][1] = $id[$i];
		}
		foreach my $i ($file_num .. $#id){
			@id_to_thread = sort { $b->[0] <=> $a->[0] } @id_to_thread if $id_to_thread[-1][0] > $id_to_thread[-2][0];
			$id_to_thread[-1][0] += $id_len{$id[$i]};
			push @{$id_to_thread[-1]}, $id[$i];
		}
	}
	foreach my $i (0 .. $#id_to_thread){
		my %seq_to_thread = ();
		foreach my $j (1 .. $#{$id_to_thread[$i]}){
			$seq_to_thread{$id_to_thread[$i][$j]} = $id_seq->{$id_to_thread[$i][$j]};
		}
		my $count = $i + 1;
		my $tmp_file_name = $fa_base_name.'.tmp'.$count;
		my $full_name = File::Spec->catfile($tmp_dir_name1, $tmp_file_name);
		store \%seq_to_thread, $full_name;
	}
	return $tmp_dir_name1;
}

#function : read search file, and store the motif length and corresponding minimum number of repeat into a hash, and return a hash reference.
sub deal_search_file{
	my $file_name = shift;
	open IN_SEARCH,'<',$file_name or die "can't open $file_name:$!";
    $_ = <IN_SEARCH>;
	die "the first line of the search file may not meet the requirement, please check the search file: $file_name\n" unless /\Am/i;
	my %molen_minum;
	$molen_minum{(split)[0]} = (split)[1] foreach <IN_SEARCH>;
    close IN_SEARCH;
	foreach my $molen (sort keys %molen_minum){
		delete $molen_minum{$molen} unless $molen_minum{$molen};
	}
	die "the contents of the search file may not meet the requirement, please check the search file: $file_name\n" if keys %molen_minum == 0;
	return \%molen_minum;
}

#function : return a hash reference that contains all SSRs of information found.
sub find_SSR_and_info{
	my ($molen_minum, $id_seq) = @_;
	my $pattern = '(';
	foreach my $molen (sort {$a <=> $b} keys %$molen_minum){
		my $minum = $molen_minum->{$molen} - 1;
		$pattern .= "([ATCG]{$molen})\\g{-1}{$minum,}|";
	}
	substr($pattern, -1) = ')';
	my %id_info;
	foreach my $id (sort keys %$id_seq){
		#this method is suitable for short sequence, so I cut each sequence into several short sequences of 1Mb in length
		my $cut_len = 1_000_000; #the length of short sequence is 1Mb
		my %shseq_loci = ();
		my $start_loci = 0;
		foreach (0 .. int( length($id_seq->{$id}) / $cut_len )){
			my $short_seq = substr $id_seq->{$id}, $start_loci, $cut_len + 10_000; #each short sequence is supplemented by another 10,000 bp, in order to prevent SSR from being cut
			$shseq_loci{$short_seq} = $start_loci + 1;
			$start_loci += $cut_len;
		}
		delete $id_seq->{$id};
		my $supple_start =  $cut_len + 1;
		my %end_loci = (); #record the end of SSRs that have been found
		foreach my $short_seq (sort { $shseq_loci{$a} <=> $shseq_loci{$b} } keys %shseq_loci){
			my $seq_start = $shseq_loci{$short_seq};
			my $front_len = $seq_start - 1;
			while(1){
				my ($SSR, @motifs) = $short_seq =~ /$pattern/;
				last unless defined $SSR;
				my ($motif) = grep defined, @motifs;
				my $SSR_len = length $SSR;
				my $motif_len = length $motif;
				my $rep_num = $SSR_len / $motif_len;
            	my $start = $-[0] + 1;
            	my $end = $+[0] + $front_len;
				substr($short_seq, 0, $start) = '';
				$start += $front_len;
				last if $start - $seq_start + 1 >= $supple_start;
				unless($end_loci{$end}){
					push @{$id_info{$id}}, [$id, $motif, $motif_len, $rep_num, $SSR_len, $start, $end];
					$end_loci{$end} = 1;
				}
				$front_len = $start;
			}
		}
	}
	return \%id_info;
	#notice : Although we have taken remedial measures, there is still a case that a very long
	####### : SSR(30,000 bp) being cut off, but this kind of situation is quite rare so we no
	####### : need to worry about it.
}

#function : use multithreads and call function find_SSR_and_info to find SSRs and their information, and return a hash reference.
sub multithreads_find_SSR_and_info{
	my ($molen_minum, $tmp_dir_name1) = @_;
	chdir $tmp_dir_name1;
	foreach my $file (glob '*'){
		my $id_seq = retrieve $file;
		my $thr = threads->create(\&find_SSR_and_info, $molen_minum, $id_seq);
	}
	my %all_SSRs;
	while(threads->list()){
    	foreach my $thr (threads->list(threads::joinable)){
        	my $SSRs = $thr->join();
			%all_SSRs = (%all_SSRs, %$SSRs);
    	}
	}
	unlink glob '*';
	chdir '..';
	rmdir $tmp_dir_name1;
	return \%all_SSRs;
}

#function : output a statistical file that contains all SSRs found to the current working directory.
sub generate_all_SSRs_stat_file{
	my ($fa_base_name, $all_SSRs) = @_;
	my $out_file_name = $fa_base_name.'.all_SSRs';
	open OUT,'>',$out_file_name or die "can't generate $out_file_name:$!";
	print OUT "id\tmotif\tmotif_len\trepeat_num\tsize\tstart\tend\n";
	foreach my $id (sort keys %$all_SSRs){
		foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$all_SSRs->{$id}}){
			my $out_row = join "\t", @$SSR;
			print OUT "$out_row\n";
		}
	}
	close OUT;
}

#----------------------------------------------------------------------------------------------------#

#function : split hash all_SSRs and id_seq into smaller hashes, respectively, and put these smaller hashes in temporary files, and return the directory name created.
sub split_SSR_and_id_seq_to_diff_file{
	my ($all_SSRs, $id_seq, $threads_num, $fa_base_name) = @_;
	my $tmp_dir_name2 = $fa_base_name.'.find_SSR_tmp2';
	if(-e $tmp_dir_name2){
		chdir $tmp_dir_name2;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name2;
	}
	mkdir $tmp_dir_name2, 0755 or die "can't create $tmp_dir_name2 directory:$!";
	my $split_file_num = ($threads_num > keys %$all_SSRs) ? keys %$all_SSRs : $threads_num;
	my $file_key_num = (keys(%$all_SSRs) % $split_file_num) ? int(keys(%$all_SSRs) / $split_file_num)+1 : keys(%$all_SSRs) / $split_file_num;
	my %out_SSRs = ();
	my %out_seqs = ();
	my ($count, $num) = (1, 1);
	foreach my $id (sort keys %$all_SSRs){
		$out_SSRs{$id} = $all_SSRs->{$id};
		$out_seqs{$id} = $id_seq->{$id};
		if(keys %out_SSRs == $file_key_num or $num == keys %$all_SSRs){
			my $new_name = $fa_base_name.'.tmp'.$count;
			my $out_file_name = File::Spec->catfile($tmp_dir_name2,$new_name);
			store [\%out_SSRs, \%out_seqs], $out_file_name;
			%out_SSRs = ();
			%out_seqs = ();
			$count++;
		}
		$num++;
	}
	return $tmp_dir_name2;
}

#function : return a hash reference that contains all SSRs of information(contains flanking sequences) found.
sub add_flanking_seq{
	my ($all_SSRs, $id_seq, $fs_len) = @_;
	my %id_SSR_info = %$all_SSRs;
	my ($left_fs, $left_fs_truelen, $right_fs, $right_fs_truelen);
	foreach my $id (sort keys %id_SSR_info){
		foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$id_SSR_info{$id}}){
			$left_fs  = ($SSR->[5] > $fs_len) ? substr $id_seq->{$id}, $SSR->[5]-$fs_len-1, $fs_len : substr $id_seq->{$id}, 0, $SSR->[5]-1;
			$left_fs_truelen  = length $left_fs;
			$right_fs = substr $id_seq->{$id}, $SSR->[6], $fs_len;
			$right_fs_truelen = length $right_fs;
			push @$SSR, ($left_fs, $left_fs_truelen, $right_fs, $right_fs_truelen);
		}
	}
	return \%id_SSR_info;
}

#function : use multithreads and call function add_flanking_seq to add flanking sequences to all SSRs found, and return a hash reference.
sub multithreads_add_flanking_seq{
	my ($fs_len, $tmp_dir_name2) = @_;
	chdir $tmp_dir_name2;
	foreach my $file (glob '*'){
		my $SSRs_and_seqs = retrieve $file;
		my $thr = threads->create(\&add_flanking_seq, $SSRs_and_seqs->[0], $SSRs_and_seqs->[1], $fs_len);
	}
	my %all_SSRs_add_fs;
	while(threads->list()){
    	foreach my $thr (threads->list(threads::joinable)){
        	my $SSRs = $thr->join();
			%all_SSRs_add_fs = (%all_SSRs_add_fs, %$SSRs);
    	}
	}
	unlink glob '*';
	chdir '..';
	rmdir $tmp_dir_name2;
	return \%all_SSRs_add_fs;
}

#function : output a statistical file that contains SSRs(add flanking sequences) information to the current working directory.
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
}

#----------------------------------------------------------------------------------------------------#

#function : return a hash reference that key is the flanking sequence, value is the frequency of existence.
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
}

#function : return a hash reference that contains the unique SSRs of flanking sequences and their information.(compared with flanking sequences of all SSRs)
sub keep_fs_uniq_SSRs{
	my ($fs_count, $all_SSRs_add_fs) = @_;
	my %fs_uni_SSRs;
	foreach my $id (sort keys %$all_SSRs_add_fs){
        foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$all_SSRs_add_fs->{$id}}){
			my $fs = $SSR->[7].'-'.$SSR->[9];
			$fs_uni_SSRs{$fs} = [ @{$SSR}[0 .. 6] ] if $fs_count->{$fs} == 1;
        }
    }
	return \%fs_uni_SSRs;
}

#function : return three hash references that contain conservative SSRs in each file and conservative SSRs statistical table, respectively.
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
	return \%fs_con_SSRs1, \%fs_con_SSRs2, \%con_compared;
}

#function : output a statistical file that contains conservative or unique SSRs to the current working directory.
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
}

#----------------------------------------------------------------------------------------------------#

#function : the values of hash id_seq are connected together and stored in a scalar, and this scalar are divided equally into smaller scalars, finally, put these smaller scalars in temporary files, and return the directory name created.
sub split_id_seq_values_to_diff_file{
	my ($id_seq, $threads_num, $fa_base_name) = @_;
	my $tmp_dir_name3 = $fa_base_name.'.find_SSR_tmp3';
	if(-e $tmp_dir_name3){
		chdir $tmp_dir_name3;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name3;
	}
	mkdir $tmp_dir_name3, 0755 or die "can't create $tmp_dir_name3 directory:$!";
    my $seqs = join '---',sort values %$id_seq;
	my $shseq_len =  int((length $seqs) / $threads_num)+1;
	my ($short_seq, $count) = ('',1);
    while($seqs){
	    $short_seq = substr $seqs, 0, $shseq_len; #notice :SSR at the boundary may be cut off, but the impact is very small
        my $new_name = $fa_base_name.'.tmp'.$count;
        my $out_file_name = File::Spec->catfile($tmp_dir_name3,$new_name);
		store \$short_seq, $out_file_name;
        $count++;
	    substr($seqs, 0, $shseq_len) = '';
    }
	return $tmp_dir_name3;
}

#function : return a array reference that contains match count of each flanking sequence compared with genome sequence.
sub get_fs_match_count{
    my ($short_seq, $new_fs_set, $fs_len) = @_;
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
	return \@match_count;
	#notice : By doing this, the flank sequence(140bp) that is not enough to set the value(150bp)
	####### : will default to mismatch. Obviously, this kind of situation rarely happen.
	#######
	####### : However, we have taken remedial measure: if the final match count is less than or equal to 2,
	####### : the flanking sequence will be regarded as unique.
}

#function : use multithreads and call function get_fs_match_count to keep the unique SSRs(the flanking sequence is unique) in the genome, and return a hash reference.
sub multithreads_keep_unique_SSRs{
	my ($tmp_dir_name3, $fs_con_SSRs, $fs_len) = @_;
	my $fs_set = [ sort {$fs_con_SSRs->{$a}[0] cmp $fs_con_SSRs->{$b}[0] or $fs_con_SSRs->{$a}[5] <=> $fs_con_SSRs->{$b}[5]} keys %$fs_con_SSRs ];
	my $new_fs_set;
	foreach my $fs (@$fs_set){
		my @fs = split /-/, $fs;
		push @{$new_fs_set}, [@fs];
	}
	unless(defined $new_fs_set){
		chdir $tmp_dir_name3;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name3;
		goto LOOP;
	}
	undef $fs_set;
	chdir $tmp_dir_name3;
	foreach my $file (glob '*'){
		my $short_seq = retrieve $file;
		my $thr = threads->create(\&get_fs_match_count, $short_seq, $new_fs_set, $fs_len);
	}
	my @all_count = ();
	while(threads->list()){
    	foreach my $thr (threads->list(threads::joinable)){
        	my $match_count = $thr->join();
			push @all_count, $match_count;
    	}
	}
	unlink glob '*';
	chdir '..';
	rmdir $tmp_dir_name3;
	my @match_count = ();
	foreach my $i (0 .. $#{$all_count[0]}){
		my $count = 0;
		$count += $all_count[$_][$i] foreach 0 .. $#all_count;
		push @match_count, $count;
	}
    my %fs_uni_SSRs;
	my $sub = 0;
	foreach my $fs (sort {$fs_con_SSRs->{$a}[0] cmp $fs_con_SSRs->{$b}[0] or $fs_con_SSRs->{$a}[5] <=> $fs_con_SSRs->{$b}[5]} keys %$fs_con_SSRs){
		$fs_uni_SSRs{$fs} = $fs_con_SSRs->{$fs} if $match_count[$sub] <= 2; #SSR with flanking sequence matching count of 0, 1, 2 is retained
		$sub++;
	}
	LOOP:
	return \%fs_uni_SSRs;
}

#----------------------------------------------------------------------------------------------------#

#function : output a statistical table that contains polymorphic SSRs candidate to the current working directory.
sub final_comparison_table{
	my ($fs_uni_SSRs1, $fs_uni_SSRs2, $con_compared, $fa_base_name1, $fa_base_name2) = @_;
	my $out_file_name = $fa_base_name1.'-and-'.$fa_base_name2.'.compare';
	open OUT,'>',$out_file_name or die "can't generate $out_file_name:$!";
	print OUT "file1_id\tfile1_motif\tfile1_repeat_num\tfile1_start\tfile1_end\tfile2_id\tfile2_motif\tfile2_repeat_num\tfile2_start\tfile2_end\tleft_fs\tleft_fs_len\tright_fs\tright_fs_len\n";
	foreach my $fs (sort {$con_compared->{$a}[0] cmp $con_compared->{$b}[0] or $con_compared->{$a}[5] <=> $con_compared->{$b}[5]} keys %$con_compared){
		if(exists $fs_uni_SSRs1->{$fs} and exists $fs_uni_SSRs2->{$fs}){
			my ($left_fs, $right_fs) = split /-/,$fs;
			my $left_fs_len  = length $left_fs;
			my $right_fs_len = length $right_fs;
			my $out_row = join "\t", (@{$con_compared->{$fs}}[0,1,3,5,6, 7,8,10,12,13], $left_fs, $left_fs_len, $right_fs, $right_fs_len);
			print OUT "$out_row\n";
		}
	}
	close OUT;
}
