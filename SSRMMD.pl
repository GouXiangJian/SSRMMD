#!/usr/bin/env perl

#date    : 2020-02-19
#author1 : Xiang-jian Gou (xjgou@stu.sicau.edu.cn)
#author2 : Hao-ran Shi (542561234@qq.com)

#load modules
use strict;
use warnings;
use 5.010;
use Cwd qw/getcwd/;
use File::Basename qw/basename/;
use File::Spec;
use Storable qw/store retrieve/;
use threads;
use Getopt::Long;

#clear buffer
$| = 1;

#get current working directory
my $cwd = getcwd;

#record version information
my $VERSION = 'SSRMMD v1.0';

#set default options
my $fasta1_name = '';
my $fasta2_name = '';
my $outDirName = 'SSRMMDOUT';
my $poly = 0;
my $motifs = '1=10,2=7,3=6,4=5,5=4,6=4';
my $minSsrLen = 10;
my $maxSsrLen = 1000;
my $fs_len = 100;
my $threads_num = 1;
my $conserMap = 'NO'; #NO NW LD
my $fsReduceFold = 0.05;
my $mismatch = 0; #0, 1 or 2
my $distanceThreshold = 0.05;
my $identityThreshold = 0.95;
my $mapScore = '1,-1,-2';
my $genome1_name = '';
my $genome2_name = '';
my $stat_file = 0;
my $inter_file = 0;
my $runLog;
my $version;
my $help;

#get options from command line
GetOptions(
		'fasta1|f1=s'    =>   \$fasta1_name,
		'fasta2|f2=s'    =>   \$fasta2_name,
		'outDir=s'       =>   \$outDirName,
		'poly=i'         =>   \$poly,
		'motifs=s'       =>   \$motifs,
		'minLen|n=i'     =>   \$minSsrLen,
		'maxLen|x=i'     =>   \$maxSsrLen,
		'length=i'       =>   \$fs_len,
		'thread=i'       =>   \$threads_num,
		'method=s'       =>   \$conserMap,
		'reduce=f'       =>   \$fsReduceFold,
		'unmatch=i'      =>   \$mismatch,
		'distance=f'     =>   \$distanceThreshold,
		'identity=f'     =>   \$identityThreshold,
		'score=s'        =>   \$mapScore,
		'genome1|g1=s'   =>   \$genome1_name,
		'genome2|g2=s'   =>   \$genome2_name,
		'stat=i'         =>   \$stat_file,
		'all=i'          =>   \$inter_file,
		'workLog=s'      =>   \$runLog,
		'version+'       =>   \$version,
		'help+'          =>   \$help,
);

#describe program information
my $usage = <<__GUIDE__;
###############################################################################
Name:
  SSRMMD - Simple Sequence Repeat Molecular Marker Developer

Author:
  Xiang-jian Gou (xjgou\@stu.sicau.edu.cn)
  Hao-ran Shi (542561234\@qq.com)

Function:
  Mining SSR loci and candidate polymorphic SSR molecular markers

Usage:
  perl SSRMMD.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  -f1 | -fasta1 <STR> : a FASTA file for mining SSR loci. (must be provided)

  -f2 | -fasta2 <STR> : another FASTA file when plan to mine polymorphic loci.

  -o  | -outDir <STR> : specify a directory for storing output file, create
                        directory if it doesn't exist. (default: SSRMMDOUT)

  -p  | -poly   <INT> : specify how to work of SSRMMD. (default: 0)
                            0 : just mining SSR loci
                            1 : mining polymorphic SSR loci

  -mo | -motifs <STR> : threshold of motif. (default: 1=10,2=7,3=6,4=5,5=4,6=4)
                           left  of equal : length of motif
                           right of equal : the minimum number of repeat

  -n  | -minLen <INT> : the minimum length of SSR. (default: 10)

  -x  | -maxLen <INT> : the maximum length of SSR. (default: 1000)

  -l  | -length <INT> : length of SSR flanking sequence. (default: 100)
                        Note: if option '-p' = 1, flanking sequences will be 
                        used to check for conservativeness and uniqueness.

  -t  | -thread <INT> : the number of threads for running. (default: 1)

  -me | -method <STR> : Algorithm for exactly checking flanking sequences
                        conservativeness. (default: NO)
                            NO : just simple check by HASH
                            LD : global alignment by Levenshtein Distance
                            NW : global alignment by Needleman-Wunsch algorithm
                        Note: 'NO' means that flanking sequences are perfectly
                        conservative, while, 'LD' and 'NW' allow mismatch or
                        indel in flanking sequences.

  -r  | -reduce   <FLOAT> : conservativeness pre-alignment by using X% flanking
                            sequences near SSR. (default: 0.05 [X% = 5%])
                            Note: -r only make sense when option '-me' = LD/NW.

  -u  | -unmatch  <INT>   : set the maximum number of mismatch base during
                            pre-alignment. (default: 0)
                                0 : no  mismatch
                                1 : one mismatch
                                2 : two mismatch
                            Note: we assume that flanking sequences near SSR
                            are highly conservative.

  -d  | -distance <FLOAT> : if option '-me' = LD, set threshold of Levenshtein
                            Distance. (default: 0.05 [5%])
                            Note: the smaller the Levenshtein Distance, the 
                            higher the sequence identity.

  -i  | -identity <FLOAT> : if option '-me' = NW, set threshold of sequence
                            identity calculated by Needleman-Wunsch algori-
                            thm. (default: 0.95 [95%])

  -sc | -score    <STR>   : mapping score of NW algorithm. (default: 1,-1,-2)
                            Note: here, 1 = match, -1 = mismatch, -2 = indel

  -g1 | -genome1  <STR>   : genome1 file of fasta1 for check flanking sequences
                            uniqueness in genome-scale. (default: fasta1 file)

  -g2 | -genome2  <STR>   : genome2 file of fasta2 for check flanking sequences
                            uniqueness in genome-scale. (default: fasta2 file)

  -st | -stat     <INT>   : whether to output SSR statistics file. (default: 0)
                                0 : not output
                                1 : yes output

  -a  | -all      <INT>   : whether to output intermediate file. (default: 0)
                                0 : not output
                                1 : yes output (be used to debug)

  -w  | -workLog  <STR>   : create a file to record run log. (default: STDOUT)

  -v  | -version          : show the version information.

  -h  | -help             : show the help information.
###############################################################################

__GUIDE__

#get map score
my ($match, $unmatch, $space) = $mapScore =~ /\A(-?\d+),(-?\d+),(-?\d+)\z/;

#check the options
die "$VERSION\n" if $version;
die $usage if $help;
die "Error: option '-f1' must be provided !\n" if ! $fasta1_name;
die "Error: option '-f2' must be provided !\n" if $poly == 1 and ! $fasta2_name;
die "Error: option '-p'  just be set to 0 or 1 !\n" if $poly != 0 and $poly != 1;
die "Error: option '-n'  must be an integer bigger than 0 !\n" if $minSsrLen <= 0 or $minSsrLen =~ /\.\d*[1-9]+/;
die "Error: option '-x'  must be an integer bigger than 0 !\n" if $maxSsrLen <= 0 or $maxSsrLen =~ /\.\d*[1-9]+/;
die "Error: option '-l'  must be an integer bigger than 0 !\n" if $fs_len <= 0 or $fs_len =~ /\.\d*[1-9]+/;
die "Error: option '-t'  must be an integer bigger than 0 !\n" if $threads_num <= 0 or $threads_num =~ /\.\d*[1-9]+/;
die "Error: option '-me' just be set to 'NO' or 'LD' or 'NW' !\n" if $conserMap ne 'NO' and $conserMap ne 'LD' and $conserMap ne 'NW';
die "Error: option '-r'  just be set to 0-0.5 !\n" if $fsReduceFold <= 0 or $fsReduceFold > 0.5;
die "Error: option '-u'  just be set to 0, 1 or 2 !\n" if $mismatch != 0 and $mismatch != 1 and $mismatch != 2;
die "Error: option '-d'  just be set to 0-0.5 !\n" if $distanceThreshold < 0 or $distanceThreshold > 0.5;
die "Error: if option '-d' be set to 0, the author suggests you that set option '-me' to 'NO' !\n" if $distanceThreshold == 0 and $conserMap ne 'NO';
die "Error: option '-i'  just be set to 0.5-1 !\n" if $identityThreshold < 0.5 or $identityThreshold > 1;
die "Error: if option '-i' be set to 1, the author suggests you that set option '-me' to 'NO' !\n" if $identityThreshold == 1 and $conserMap ne 'NO';
die "Error: option '-sc' have a problem !\n" if (! defined($match)) or (! defined($unmatch)) or (! defined($space));
die "Error: option '-st' just be set to 0 or 1 !\n" if $stat_file != 0 and $stat_file != 1;
die "Error: option '-a'  just be set to 0 or 1 !\n" if $inter_file != 0 and $inter_file != 1;

#specify conservativeness threshold
my $conserThreshold = $conserMap eq 'LD' ? $distanceThreshold : $identityThreshold;

#specify genome file
$genome1_name = $fasta1_name if ! $genome1_name;
$genome2_name = $fasta2_name if $fasta2_name and ! $genome2_name;

#get base name of fasta file
my $fa_base_name1 = basename $fasta1_name;
my $fa_base_name2 = basename $fasta2_name if $fasta2_name;

#create output directory
mkdir $outDirName, 0755 or die "Error: can't create directory '$outDirName' : $!" if ! -d $outDirName;
chdir $outDirName or die "Error: can't cd to directory '$outDirName' : $!";

#create log file
my $log;
open $log, '>', $runLog or die "Error: can't generate $runLog : $!" if defined $runLog;
my $logMark = defined $runLog ? $runLog : 'STDOUT';
my $optionSet = <<__OPTION__;
parameter settings:
  -f1 $fasta1_name
  -f2 $fasta2_name
  -o  $outDirName
  -p  $poly
  -mo $motifs
  -n  $minSsrLen
  -x  $maxSsrLen
  -l  $fs_len
  -t  $threads_num
  -me $conserMap
  -r  $fsReduceFold
  -u  $mismatch
  -d  $distanceThreshold
  -i  $identityThreshold
  -sc $mapScore
  -g1 $genome1_name
  -g2 $genome2_name
  -st $stat_file
  -a  $inter_file
  -w  $logMark

run record:
__OPTION__

#write options information to logfile
writeLogInfo($log, $optionSet);



####################################################################################################
#                                      start main program ...                                      #
####################################################################################################



#Module1 : mining SSR loci

if ($poly == 0) {

	#record run step

	my $localTime1 = localtime;

	writeLogInfo($log, "  Start to mine SSR loci for fasta1 : $localTime1\n");


	#get some base information

	my $fa_base_name = basename $fasta1_name;

	my $molen_minum = get_motifs_info($motifs);

	my $id_seq = deal_fa_file($fasta1_name, $outDirName, $cwd);


	#split seq to tmp file

	my ($tmp_dir_name1, $final_threads_num, $id_len) = split_id_seq_to_diff_file($id_seq, $threads_num, $fa_base_name);

	undef $id_seq;


	#mining SSR loci

	my $all_SSRs = multithreads_find_SSR_and_info($molen_minum, $tmp_dir_name1, $final_threads_num, $fs_len, $minSsrLen, $maxSsrLen);

	undef $tmp_dir_name1;


	#output statistics information

	some_simple_statistics($fa_base_name, $all_SSRs, $id_len, $molen_minum) if $stat_file;

	undef $id_len;


	#output SSR information

	generate_SSRs_with_fs_stat_file($fa_base_name, '.SSRs', $all_SSRs);

	undef $all_SSRs;


	#record run step

	my $localTime2 = localtime;

	writeLogInfo($log, "  Run finish : $localTime2\n");
}

#Module2 : developing candidate polymorphic SSR

else {

	#declare an array for storing some tmp information

	my @retrieve_files;

	my $countNum = 1;


	#================================================================================

	#step1 : mining SSRs loci


	foreach my $fa_file_name ($fasta1_name, $fasta2_name) {

		#record run step

		my $localTime1 = localtime;

		writeLogInfo($log, "  Start to mine SSR loci for fasta$countNum : $localTime1\n");

		$countNum++;


		#get some base information

		my $fa_base_name = basename $fa_file_name;

		my $molen_minum = get_motifs_info($motifs);

		my $id_seq = deal_fa_file($fa_file_name, $outDirName, $cwd);


		#split seq to tmp file

		my ($tmp_dir_name1, $final_threads_num, $id_len) = split_id_seq_to_diff_file($id_seq, $threads_num, $fa_base_name);

		push @retrieve_files, $id_len;

		undef $id_seq;


		#mining SSR loci

		my $all_SSRs = multithreads_find_SSR_and_info($molen_minum, $tmp_dir_name1, $final_threads_num, $fs_len, $minSsrLen, $maxSsrLen);

		undef $tmp_dir_name1;


		#output SSR information

		generate_SSRs_with_fs_stat_file($fa_base_name, '.SSRs', $all_SSRs);


		#output statistics information

		some_simple_statistics($fa_base_name, $all_SSRs, $id_len, $molen_minum) if $stat_file;

		undef $id_len;


		#store SSR information to tmp file

		push @retrieve_files, StoreToFile($fa_base_name, 'ssr', $all_SSRs);

		undef $all_SSRs;
	}

		#================================================================================

		#step2 : keep all flanking sequence conservative SSRs


		#record run step

		my $localTime2 = localtime;

		writeLogInfo($log, "  Start to check conservativeness of flanking sequences : $localTime2\n");


		#deal fasta1 to exclude non-unique SSR by myself

		my $all_SSRs_add_fs1 = RetrieveFromFile($retrieve_files[1]);

		my $fs_count1 = generate_fs_list($all_SSRs_add_fs1);

		my $fs_uni_SSRs1 = keep_fs_uniq_SSRs($fs_count1, $all_SSRs_add_fs1);

		undef $fs_count1;

		undef $all_SSRs_add_fs1;


		#deal fasta2 to exclude non-unique SSR by myself

		my $all_SSRs_add_fs2 = RetrieveFromFile($retrieve_files[3]);

		my $fs_count2 = generate_fs_list($all_SSRs_add_fs2);

		my $fs_uni_SSRs2 = keep_fs_uniq_SSRs($fs_count2, $all_SSRs_add_fs2);

		undef $fs_count2;

		undef $all_SSRs_add_fs2;


		#simple conservativeness check by HASH

		my ($fs_con_SSRs1, $fs_con_SSRs2, $con_compared) = keep_fs_cons_SSRs($fs_uni_SSRs1, $fs_uni_SSRs2);


		#exactly conservativeness check by LD or NW

		if ($conserMap ne 'NO') {

			#1 store conservativeness information of compare to tmp file

			my $map_cmp_file = StoreToFile($fa_base_name1, 'mapping-cmp', $con_compared);

			undef $con_compared;


			#2 get rest pre-unique SSRs of fasta1 from all unique sets1

			my $restUniSSRs1 = getRestSSRsFromUniSet($fs_uni_SSRs1, $fs_con_SSRs1, $fsReduceFold);

			undef $fs_uni_SSRs1;

			my $map_con_file1 = StoreToFile($fa_base_name1, 'mapping', $fs_con_SSRs1);

			undef $fs_con_SSRs1;


			#3 get rest pre-unique SSRs of fasta2 from all unique sets2

			my $restUniSSRs2 = getRestSSRsFromUniSet($fs_uni_SSRs2, $fs_con_SSRs2, $fsReduceFold);

			undef $fs_uni_SSRs2;

			my $map_con_file2 = StoreToFile($fa_base_name2, 'mapping', $fs_con_SSRs2);

			undef $fs_con_SSRs2;


			#4 pre-alignment for checking conservativeness by short sequence

			my $common = prealign_highly_conservativeness($restUniSSRs1, $restUniSSRs2, $mismatch);

			undef $restUniSSRs1;

			undef $restUniSSRs2;


			#5 split information to tmp file

			my $base_name = basename $fasta1_name;

			my $tmp_dir_name2 = splitRestSSRsToTempFiles($common, $threads_num, $base_name);

			undef $common;


			#6 start to conservativeness mapping by LD or NW

			my ($con_SSRs1, $con_SSRs2, $compare) = multiThreadsMapFsToFindConSSR($tmp_dir_name2, $conserThreshold, $conserMap);


			#7 get conservativeness fasta1 file

			$fs_con_SSRs1 = RetrieveFromFile($map_con_file1);

			$fs_con_SSRs1 = mergeScalar($fs_con_SSRs1, $con_SSRs1);
		
			undef $con_SSRs1;


			#8 get conservativeness fasta2 file

			$fs_con_SSRs2 = RetrieveFromFile($map_con_file2);

			$fs_con_SSRs2 = mergeScalar($fs_con_SSRs2, $con_SSRs2);

			undef $con_SSRs2;


			#9 get conservativeness compare file

			$con_compared = RetrieveFromFile($map_cmp_file);

			$con_compared = mergeScalar($con_compared, $compare);

			undef $compare;
		}


		#free memory

		undef $fs_uni_SSRs1;

		undef $fs_uni_SSRs2;


		#output conservativeness file of fasta1 for debug

		generate_cons_or_uniq_SSRs_stat_file($fa_base_name1, '.conservative', $fs_con_SSRs1) if $inter_file;
		

		#output conservativeness file of fasta2 for debug

		generate_cons_or_uniq_SSRs_stat_file($fa_base_name2, '.conservative', $fs_con_SSRs2) if $inter_file;


		#store conservativeness information of fasta1 to tmp file

		my $fs_con_SSRs_file1 = StoreToFile($fa_base_name1, 'con', $fs_con_SSRs1);
		
		undef $fs_con_SSRs1;

		
		#store conservativeness information of fasta2 to tmp file

		my $fs_con_SSRs_file2 = StoreToFile($fa_base_name2, 'con', $fs_con_SSRs2);
		
		undef $fs_con_SSRs2;


		#store conservativeness information of compare to tmp file

		my $con_compare_file = StoreToFile("$fa_base_name1-$fa_base_name2", 'cmp', $con_compared);
		
		undef $con_compared;


		#================================================================================

		#step3 : keep all unique SSRs


		#record run step

		my $localTime3 = localtime;

		writeLogInfo($log, "  Start to check uniqueness of flanking sequences for fasta1 : $localTime3\n");


		#deal fasta1 to check uniqueness of flanking sequences

		my $tmp_dir_name3_1 = split_id_seq_values_to_diff_file($genome1_name, $outDirName, $cwd, $threads_num, $retrieve_files[0]);

		undef $retrieve_files[0]; 

		my $fs_uniq_SSRs1 = multithreads_keep_unique_SSRs($tmp_dir_name3_1, $fs_con_SSRs_file1, $fs_len);

		undef $tmp_dir_name3_1;

		undef $fs_con_SSRs_file1;


		#output uniqueness file of fasta1 for debug

		generate_cons_or_uniq_SSRs_stat_file($fa_base_name1, '.unique', $fs_uniq_SSRs1) if $inter_file;


		#store uniqueness information of fasta1 to tmp file

		my $uniq_SSRs_file1 = StoreToFile($fa_base_name1, 'uni', $fs_uniq_SSRs1);

		undef $fs_uniq_SSRs1;


		#record run step

		my $localTime4 = localtime;

		writeLogInfo($log, "  Start to check uniqueness of flanking sequences for fasta2 : $localTime4\n");


		#deal fasta2 to check uniqueness of flanking sequences

		my $tmp_dir_name3_2 = split_id_seq_values_to_diff_file($genome2_name, $outDirName, $cwd, $threads_num, $retrieve_files[2]);

		undef $retrieve_files[2];

		my $fs_uniq_SSRs2 = multithreads_keep_unique_SSRs($tmp_dir_name3_2, $fs_con_SSRs_file2, $fs_len);
		
		undef $tmp_dir_name3_2;

		undef $fs_con_SSRs_file2;


		#output uniqueness file of fasta2 for debug

		generate_cons_or_uniq_SSRs_stat_file($fa_base_name2, '.unique', $fs_uniq_SSRs2) if $inter_file;


		#retrieve uniqueness information of fasta1

		$fs_uniq_SSRs1 = RetrieveFromFile($uniq_SSRs_file1);


		#retrieve uniqueness information of compare

		$con_compared = RetrieveFromFile($con_compare_file);


		#================================================================================

		#step4 : output a table of candidate polymorphic SSRs


		#record run step

		my $localTime5 = localtime;

		writeLogInfo($log, "  Start to do polymorphism check : $localTime5\n");


		#final compare

		final_comparison_table($fs_uniq_SSRs1, $fs_uniq_SSRs2, $con_compared, $fa_base_name1, $fa_base_name2, $conserMap);


		#record run step

		my $localTime6 = localtime;

		writeLogInfo($log, "  Run finish : $localTime6\n");
}



####################################################################################################
#                                      end main program !!!                                        #
####################################################################################################



#here are all subroutines are used by this program:

#01 writeLogInfo
#02 StoreToFile
#03 RetrieveFromFile
#04 deal_fa_file
#05 split_id_seq_to_diff_file
#06 average_allocation_algorithm
#07 get_motifs_info
#08 is_false_motif
#09 find_SSR_and_info
#10 multithreads_find_SSR_and_info
#11 some_simple_statistics
#12 generate_SSRs_with_fs_stat_file
#13 generate_fs_list
#14 keep_fs_uniq_SSRs
#15 keep_fs_cons_SSRs
#16 getRestSSRsFromUniSet
#17 prealign_highly_conservativeness
#18 mergeScalar
#19 mapFsToFindConSSR
#20 multiThreadsMapFsToFindConSSR
#21 splitRestSSRsToTempFiles
#22 Levenshtein_Distance
#23 Needleman_Wunsch
#24 generate_cons_or_uniq_SSRs_stat_file
#25 split_id_seq_values_to_diff_file
#26 get_fs_match_count
#27 sldingWindows
#28 multithreads_keep_unique_SSRs
#29 final_comparison_table

#======================================================================


#function1 : record run information
sub writeLogInfo {
	my ($handle, $info) = @_;
	if ($handle) {
		print $handle $info;
	}
	else {
		print $info;
	}
}


#function2 : store information to tmp file, and return a file name
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


#function3 : retrieve information from tmp file, and return the information
sub RetrieveFromFile {
	my $file = shift;	
	die "Error: no find file '$file' : $!" unless -e $file;
	my $info = retrieve $file;
	unlink $file;
	undef $file;
	return $info;
}


#function4 : read fasta file, and store the sequence name and corresponding sequence into a hash, and return a hash reference.
sub deal_fa_file {
	my ($file_name, $dirName, $cwd) = @_;
	chdir $cwd;
	open IN_FA, '<', $file_name or die "Error: can't open $file_name : $!";
	my ($id, %id_seq) = ('', ());
	while(<IN_FA>){
		s/[\r\n]+//;
		/\A>(.*)/ ? ( $id = $1 =~ s/\s/_/gr ) : ( $_ = uc, s/[^A-Z]//g, $id_seq{$id} .= $_ );
	}
	close IN_FA;
	die "the format of fasta file may be incorrect, please check the fasta file: $file_name\n" if keys %id_seq == 0;
	chdir $dirName;
	return \%id_seq;
}


#function5 : split hash id_seq into smaller hashes, and put these smaller hashes in temporary files.
sub split_id_seq_to_diff_file {
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
	mkdir $tmp_dir_name1, 0755 or die "Error: can't create '$tmp_dir_name1' directory:$!";
	my %id_len;
	$id_len{$_} = length $id_seq->{$_} foreach keys %$id_seq;
	my $seq_sum = scalar keys %id_len;
	my $final_threads_num = ($threads_num >= $seq_sum) ? $seq_sum : $threads_num;
	my @id_to_thread = average_allocation_algorithm($final_threads_num, %id_len);
	if ($seq_sum == 1) {
		my $full_name = File::Spec->catfile($tmp_dir_name1, "1-1");
		store $id_seq, $full_name;
	}
	else {
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
				#delete $id_seq->{$id_to_thread[$i][$j]};
				$each_file_size += $id_len{$id_to_thread[$i][$j]};
			}
			my $full_name = File::Spec->catfile($tmp_dir_name1, "$count-$in_thread_count");
			store \%seq_to_thread, $full_name;
			undef %seq_to_thread;
		}
	}
	undef $id_seq;
	undef @id_to_thread;
	return $tmp_dir_name1, $final_threads_num, \%id_len;
}


#function6 : allocate the sequence evenly to different threads by length, and return a two-dimensional array
sub average_allocation_algorithm {
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
	undef @id;
	undef %id_len;
	return @id_to_thread;
}


#function7 : store the motif length and corresponding minimum number of repeat into a hash, and return a hash reference.
sub get_motifs_info {
	my $motifs_info = shift;
	my @digit = $motifs_info =~ /(\d+)=(\d+)/g;
	die "Error: option '-mo' have a problem !\n" if @digit == 0 or @digit % 2; 
	my %molen_minum = @digit;
	return \%molen_minum;
}


#function8 : judge if motif is false, 1 = false, 0 = true
sub is_false_motif {
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
            my $sMotif = substr $tmp, 0, $len;
            $motif_sub{$sMotif} = 1;
            substr($tmp, 0, $len) = '';
        }
        if(keys %motif_sub == 1) {
            $judge = 1;
            last;
        }
    }
    return $judge;
}


#function9 : return a hash reference that contains all SSRs of information found.
sub find_SSR_and_info {
	my ($molen_minum, $files, $fs_len, $minSsrLen, $maxSsrLen) = @_;
	my %id_info;
	foreach my $file (@$files) {
		my $id_seq = retrieve $file;
		foreach my $id (sort keys %$id_seq) {
			foreach my $molen (sort {$a <=> $b} keys %$molen_minum) {
				my $remain = $molen_minum->{$molen} - 1;
				while ($id_seq->{$id} =~ /(([ATCG]{$molen})\g{-1}{$remain,})/g) {
					my ($SSR, $motif) = ($1, $2);
					my $SSR_len = length $SSR;
					my $rep_num = $SSR_len / $molen;
					my $start = $-[0] + 1;
					my $end   = $+[0];
					pos($id_seq->{$id}) += 1 - $molen; #backtrack (makesure adequate mining of composite SSR)
					if (! is_false_motif($motif)) {
						if ($SSR_len >= $minSsrLen and $SSR_len <= $maxSsrLen) {
							my $left_fs = '';
							my $left_fs_truelen = 0;
							my $right_fs = '';
							my $right_fs_truelen = 0;
							$left_fs  = ($start > $fs_len) ? substr $id_seq->{$id}, $start-$fs_len-1, $fs_len : substr $id_seq->{$id}, 0, $start-1;
							$left_fs_truelen  = length $left_fs;
							$right_fs = substr $id_seq->{$id}, $end, $fs_len;
							$right_fs_truelen = length $right_fs;
							push @{$id_info{$id}}, [$id, $motif, $molen, $rep_num, $SSR_len, $start, $end, $left_fs, $left_fs_truelen, $right_fs, $right_fs_truelen];
						}
					}
				}
				pos($id_seq->{$id}) = 0; #reset the position
			}
			delete $id_seq->{$id};
		}
		undef $id_seq;
	}
	return \%id_info;
}


#function10 : use multithreads and call function find_SSR_and_info to find SSRs and their information, and return a hash reference.
sub multithreads_find_SSR_and_info {
	my ($molen_minum, $tmp_dir_name1, $final_threads_num, $fs_len, $minSsrLen, $maxSsrLen) = @_;
	my %all_SSRs;
	if ($final_threads_num == 1) {
		my @files = glob "$tmp_dir_name1/1-*";
		die "Error: can't find any file in tmp directory '$tmp_dir_name1' : $!" if @files == 0;
		my $SSRs = find_SSR_and_info($molen_minum, \@files, $fs_len, $minSsrLen, $maxSsrLen);
		%all_SSRs = %$SSRs;
		undef $SSRs;
	}
	else {
		foreach my $i (1 .. $final_threads_num) {
			my @files = glob "$tmp_dir_name1/$i-*";
			die "Error: can't find any file with prefix $i in tmp directory '$tmp_dir_name1' : $!" if @files == 0;
			my $thr = threads->create(\&find_SSR_and_info, $molen_minum, \@files, $fs_len, $minSsrLen, $maxSsrLen);
		}
		while(threads->list()){
    		foreach my $thr (threads->list(threads::joinable)){
				my $SSRs = $thr->join();
				%all_SSRs = (%all_SSRs, %$SSRs);
				undef $SSRs;
    		}
		}
	}
	if (-e $tmp_dir_name1) {
		chdir $tmp_dir_name1;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name1;
	}
	undef $tmp_dir_name1;
	return \%all_SSRs;
}


#function11 : output some simple statistics about all SSRs loci
sub some_simple_statistics {
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
	open OUT_ALL, '>', $out_file_name or die "Error: can't generate $out_file_name : $!";
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
	undef $all_SSRs;
	undef %id_len;
	undef $id_len;
	undef %id_sum;
	undef %motif_sum;
	undef %motif_num;
}


#function12 : output a statistical file that contains SSRs(add flanking sequences) information to the current working directory.
sub generate_SSRs_with_fs_stat_file {
	my ($fa_base_name, $out_file_suffix, $SSRs_info) = @_;
	my $out_file_name = $fa_base_name.$out_file_suffix;
	open my $OUT,'>',$out_file_name or die "Error: can't generate $out_file_name:$!";
	print $OUT "number\tid\tmotif\tmotif_length\trepeat_number\tsize\tstart\tend\tleft_fs\tleft_fs_length\tright_fs\tright_fs_length\n";
	my $count = 1;
	foreach my $id (sort keys %$SSRs_info){
		foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$SSRs_info->{$id}}){
			my $out_row = join "\t", ($count, @$SSR);
			print $OUT "$out_row\n";
			$count++;
		}
	}
	close $OUT;
	undef $SSRs_info;
}


#function13 : return a hash reference that key is the flanking sequence, value is the frequency of existence.
sub generate_fs_list {
	my $all_SSRs_add_fs = shift;
	my %fs_count;
	foreach my $id (sort keys %$all_SSRs_add_fs){
        foreach my $SSR (sort {$a->[5] <=> $b->[5]} @{$all_SSRs_add_fs->{$id}}){
			my $fs = $SSR->[7].'-'.$SSR->[9]; #connect the left flanking sequence and right flanking sequence together, and separate them by using '-'
			$fs_count{$fs}++;
        }
    }
	undef $all_SSRs_add_fs;
	return \%fs_count;
}


#function14 : return a hash reference that contains the unique SSRs of flanking sequences and their information.(compared with flanking sequences of all SSRs)
sub keep_fs_uniq_SSRs {
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
sub keep_fs_cons_SSRs {
	my ($fs_uni_SSRs1, $fs_uni_SSRs2) = @_;
	my (%fs_con_SSRs1, %fs_con_SSRs2, %con_compared);
	foreach my $fs1 (sort keys %$fs_uni_SSRs1){
		if(exists $fs_uni_SSRs2->{$fs1} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$fs1}[1]){ #the flanking sequence is found in every file, and SSR motif is the same
			$fs_con_SSRs1{$fs1} = $fs_uni_SSRs1->{$fs1};
			$fs_con_SSRs2{$fs1} = $fs_uni_SSRs2->{$fs1};
			$con_compared{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$fs1}}, 'same', 'NA', 'NA', 'NA', 'NA'];
		}
	}
	undef $fs_uni_SSRs1;
	undef $fs_uni_SSRs2;
	return \%fs_con_SSRs1, \%fs_con_SSRs2, \%con_compared;
}


#function16 : exclude SSR that already is conservative SSR, and return a hash reference.
sub getRestSSRsFromUniSet {
	my ($fs_uni_SSRs, $fs_con_SSRs, $reduceFactor) = @_;
	my (%restSet, %mark);
	foreach my $fs (sort keys %$fs_uni_SSRs) {
		if (exists $fs_con_SSRs->{$fs}) {
			delete $fs_uni_SSRs->{$fs};
		}
		else {
			my ($left_fs, $right_fs) = split /-/, $fs;
			next unless $left_fs and $right_fs;
			my $leftFsLen = length $left_fs;
			my $leftIntLen = int($leftFsLen*$reduceFactor);
			my $leftFsReduceLen = $leftIntLen + ($leftFsLen*$reduceFactor > $leftIntLen ? 1 : 0);
			my $leftFsReduce = substr($left_fs, -$leftFsReduceLen);
			my $rightFsLen = length $right_fs;
			my $rightIntLen = int($rightFsLen*$reduceFactor);
			my $rightFsReduceLen = $rightIntLen + ($rightFsLen*$reduceFactor > $rightIntLen ? 1 : 0);
			my $rightFsReduce =  substr($right_fs, 0, $rightFsReduceLen);
			my $short_fs = $leftFsReduce.'-'.$rightFsReduce;
			if (! exists $restSet{$short_fs}) {
				if (! exists $mark{$short_fs}) {
					push @{$fs_uni_SSRs->{$fs}}, $left_fs, $right_fs;
					$restSet{$short_fs} = $fs_uni_SSRs->{$fs};
				}
			}
			else {
				delete $restSet{$short_fs};
				$mark{$short_fs} = 1;
			}
		}
	}
	undef $fs_uni_SSRs;
	undef $fs_con_SSRs;
	undef %mark;
	return \%restSet;
}


#function17 : pre-alignment by highly conservativeness, and return a hash reference.
sub prealign_highly_conservativeness {
	my ($fs_uni_SSRs1, $fs_uni_SSRs2, $judge) = @_; # $judge just be 0, 1 or 2
	my %common = ();
	if ($judge == 0) {
		foreach my $fs1 (sort keys %$fs_uni_SSRs1) {
			if (exists $fs_uni_SSRs2->{$fs1} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$fs1}[1]) {
				$common{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$fs1}} ];
			}
		}
	}
	elsif($judge == 1) {
		foreach my $fs1 (sort keys %$fs_uni_SSRs1) {
			if (exists $fs_uni_SSRs2->{$fs1} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$fs1}[1]) {
				$common{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$fs1}} ];
			}
			else {
				my $len = length $fs1;
				THIS1:foreach my $base (qw/A G C T/) {
					foreach my $i (0 .. $len-1) {
						my $seq = $fs1;
						substr($seq, $i, 1) = $base;
						if (exists $fs_uni_SSRs2->{$seq} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$seq}[1]) {
							$common{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$seq}} ];
							last THIS1;
						}
					}
				}
			}
		}
	}
	elsif ($judge == 2) {
		foreach my $fs1 (sort keys %$fs_uni_SSRs1) {
			if (exists $fs_uni_SSRs2->{$fs1} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$fs1}[1]) {
				$common{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$fs1}} ];
			}
			else {
				my $len1 = length $fs1;
				my $mark = 0;
				THIS2:foreach my $base (qw/A G C T/) {
					foreach my $i (0 .. $len1-1) {
						my $seq = $fs1;
						substr($seq, $i, 1) = $base;
						if (exists $fs_uni_SSRs2->{$seq} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$seq}[1]) {
							$common{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$seq}} ];
							$mark = 1;
							last THIS2;
						}
					}
				}
				if ($mark == 0) {
					my $len2 = length $fs1;
					THIS3:foreach my $base1 (qw/A G C T/) {
						foreach my $base2 (qw/A G C T/) {
							foreach my $i (0 .. $len2-1) {
								foreach my $j ($i+1 .. $len2-1) {
									my $seq = $fs1;
									substr($seq, $i, 1) = $base1;
									substr($seq, $j, 1) = $base2;
									if (exists $fs_uni_SSRs2->{$seq} and $fs_uni_SSRs1->{$fs1}[1] eq $fs_uni_SSRs2->{$seq}[1]) {
										$common{$fs1} = [ @{$fs_uni_SSRs1->{$fs1}}, @{$fs_uni_SSRs2->{$seq}} ];
										last THIS3;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {
		die "Error: the program have a bug !\n";
	}
	undef $fs_uni_SSRs1;
	undef $fs_uni_SSRs2;
	return \%common;
}


#function18 : merge two scalar into a single scalar, and return a hash reference.
sub mergeScalar {
    my ($one, $two) = @_;
    foreach my $each (sort keys %$two) {
        $one->{$each} = $two->{$each};
    }
    undef $two;
    return $one;
}


#function19 : map two seqs for finding conservative SSRs, and return a array reference.
sub mapFsToFindConSSR {
	my ($file, $conserThreshold, $method) = @_;
	my $restSSRs = retrieve $file;
	undef $file;
	my $codeRef;
	if ($method eq 'LD') {
		$codeRef = \&Levenshtein_Distance;
	}
	elsif ($method eq 'NW') {
		$codeRef = \&Needleman_Wunsch;
	}
	my (%con_SSRs1, %con_SSRs2, %compare);
	foreach my $fs (sort keys %$restSSRs) {
		my $left_fs1  = $restSSRs->{$fs}[7];
		my $right_fs1 = $restSSRs->{$fs}[8];
		my $left_fs2  = $restSSRs->{$fs}[16];
		my $right_fs2 = $restSSRs->{$fs}[17];
		my ($judge1, $Lresult, $Lwhich) = $codeRef->($left_fs1, $left_fs2, $conserThreshold); # map left seq
		if ($judge1) {
            my ($judge2, $Rresult, $Rwhich) = $codeRef->($right_fs1, $right_fs2, $conserThreshold); # map right seq
            if ($judge2) {
				my $real_fs1 = $left_fs1.'-'.$right_fs1;
				$con_SSRs1{$real_fs1} = [ @{$restSSRs->{$fs}}[0 .. 6] ];
				my $real_fs2 = $left_fs2.'-'.$right_fs2;
				$con_SSRs2{$real_fs2} = [ @{$restSSRs->{$fs}}[9 .. 15] ];
                $compare{$real_fs1} = [ @{$con_SSRs1{$real_fs1}}, @{$con_SSRs2{$real_fs2}}, $real_fs2, $Lresult, $Lwhich, $Rresult, $Rwhich ];
        	}
        }
    }
    undef $restSSRs;
    my @outCont = (\%con_SSRs1, \%con_SSRs2, \%compare);
    return \@outCont;
}


#function20 : use multithreads to map two seqs for finding conservative SSRs, and return three hash references.
sub multiThreadsMapFsToFindConSSR {
	my ($tmp_dir_name, $conserThreshold, $method) = @_;
	my (%con_SSRs1, %con_SSRs2, %compare);
	my @files = glob "$tmp_dir_name/*";
	if (@files == 1) {
		my $outCont = mapFsToFindConSSR($files[0], $conserThreshold, $method);
		%con_SSRs1 = %{$outCont->[0]};
		%con_SSRs2 = %{$outCont->[1]};
		%compare   = %{$outCont->[2]};
		undef $outCont;
	}
	elsif (@files >= 2) {
		foreach my $file (@files) {
			my $thr = threads->create(\&mapFsToFindConSSR, $file, $conserThreshold, $method);
		}
		while(threads->list()){
			foreach my $thr (threads->list(threads::joinable)){
				my $outCont = $thr->join();
				%con_SSRs1 = (%con_SSRs1, %{$outCont->[0]});
				%con_SSRs2 = (%con_SSRs2, %{$outCont->[1]});
				%compare   = (%compare,   %{$outCont->[2]});
				undef $outCont;
			}
		}
	}
	else {
		die "Error: can't find any file in tmp directory '$tmp_dir_name' : $!";
	}
	if (-e $tmp_dir_name) {
		chdir $tmp_dir_name;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name;
	}
	undef $tmp_dir_name;
	return \%con_SSRs1, \%con_SSRs2, \%compare;
}


#function21 : split hash into smaller hashes, and put these smaller hashes in temporary files.
sub splitRestSSRsToTempFiles {
	my ($restSSRs, $threads_num, $fa_base_name) = @_;
	my ($sec, $min, $hour, $day, $mon, $year, undef, undef, undef) = localtime;
	$mon  += 1;
	$year += 1900;
	my $time = "$year$mon$day-$hour$min$sec";
	my $tmp_dir_name = $fa_base_name."_$time.SSRMMD_tmp2";
	if (-e $tmp_dir_name) {
		chdir $tmp_dir_name;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name;
	}
	mkdir $tmp_dir_name, 0755 or die "Error: can't create '$tmp_dir_name' directory : $!";
	if ($threads_num == 1) {
		my $full_name = File::Spec->catfile($tmp_dir_name, "thread.1");
		store $restSSRs, $full_name;
	}
	else {
		my $total = scalar keys %$restSSRs;
		my $eachThreadSsrNum = int($total/$threads_num) + ($total % $threads_num ? 1 : 0); 
		my %outThreadSSRs;
		my $fileNum = 1;
		my $count = 0;
		foreach my $each (sort keys %$restSSRs) {
			if ($count < $eachThreadSsrNum) {
				$outThreadSSRs{$each} = $restSSRs->{$each};
				$count++;
			}
			else {
				my $full_name = File::Spec->catfile($tmp_dir_name, "thread.$fileNum");
				store \%outThreadSSRs, $full_name;
				$fileNum++;
				%outThreadSSRs = ();
				$outThreadSSRs{$each} = $restSSRs->{$each};
				$count = 1;
			}
		}
		if ($count != 0) {
			my $full_name = File::Spec->catfile($tmp_dir_name, "thread.$fileNum");
			store \%outThreadSSRs, $full_name;
		}
		undef %outThreadSSRs;
	}
	undef $restSSRs;
	return $tmp_dir_name;
}


#function22 : global alignment by using levenshtein distance.
sub Levenshtein_Distance {
    my ($one, $two, $threshold) = @_;
    my $seq1 = [ split //, $one ];
    my $seq1Len = @$seq1;
    undef $one;
    my $seq2 = [ split //, $two ];
    my $seq2Len = @$seq2;
    undef $two;
    my @matrix_score;
    $matrix_score[0][0] = 0;
    foreach my $j (1 .. $seq2Len) {
        $matrix_score[0][$j] = $matrix_score[0][$j-1] + 1;
    }
    foreach my $i (1 .. $seq1Len) {
        $matrix_score[$i][0] = $matrix_score[$i-1][0] + 1;
    }
    foreach my $i (1 .. $seq1Len) {
        foreach my $j (1 .. $seq2Len) {
            my $road_1_score = $matrix_score[$i][$j-1] + 1; #reach each cell from the left
            my $road_2_score = $matrix_score[$i-1][$j] + 1; #reach each cell from the up
            my $mark = $seq1->[$i-1] eq $seq2->[$j-1] ? 0 : 1;
            my $road_3_score = $matrix_score[$i-1][$j-1] + $mark; #reach each cell from the upper left
            if ($road_1_score > $road_2_score) {
                $matrix_score[$i][$j] = $road_2_score > $road_3_score ? $road_3_score : $road_2_score;
            }
            else{
                $matrix_score[$i][$j] = $road_1_score > $road_3_score ? $road_3_score : $road_1_score;
            }
        }
    }
    my $maxLen = $seq1Len > $seq2Len ? $seq1Len : $seq2Len;
    my $disRatio = sprintf "%.3f", $matrix_score[-1][-1]/$maxLen;
	my $judge = $disRatio <= $threshold ? 1 : 0;
	undef @matrix_score;
    return $judge, $disRatio, 'D';
}


#function23 : global alignment by using Needleman-Wunsch algorithm
sub Needleman_Wunsch {
    my ($one, $two, $threshold) = @_;
    my $seq1 = [ split //, $one ];
    my $seq1Len = @$seq1;
    undef $one;
    my $seq2 = [ split //, $two ];
    my $seq2Len = @$seq2;
    undef $two;
    my (@matrix_score, @matrix_arrow); #in arrow matrix, -1 mean left, 1 mean up, and 0 mean upper left
    $matrix_score[0][0] = 0;
    $matrix_arrow[0][0] = 0;
    foreach my $j (1 .. $seq2Len) {
        $matrix_score[0][$j] = $matrix_score[0][$j-1] + $space;
        $matrix_arrow[0][$j] = -1;
    }
    foreach my $i (1 .. $seq1Len) {
        $matrix_score[$i][0] = $matrix_score[$i-1][0] + $space;
        $matrix_arrow[$i][0] = 1;
    }
    foreach my $i (1 .. $seq1Len) {
        foreach my $j (1 .. $seq2Len) {
            my $road_1_score = $matrix_score[$i][$j-1] + $space; #reach each cell from the left
            my $road_2_score = $matrix_score[$i-1][$j] + $space; #reach each cell from the up
            my $mark = $seq1->[$i-1] eq $seq2->[$j-1] ? $match : $unmatch;
            my $road_3_score = $matrix_score[$i-1][$j-1] + $mark; #reach each cell from the upper left
            if ($road_1_score > $road_2_score) {
                $matrix_score[$i][$j] = $road_1_score > $road_3_score ? $road_1_score : $road_3_score;
            }
            else{
                $matrix_score[$i][$j] = $road_2_score > $road_3_score ? $road_2_score : $road_3_score;
            }
            #if there are multiple paths for the maximum score, road 3 will be used first, road 2 will be used second, and road 1 will be used last
            $matrix_arrow[$i][$j] = $matrix_score[$i][$j] == $road_3_score ? 0 : ($matrix_score[$i][$j] == $road_2_score ? 1 : -1);
        }
    }
    #backtracking
    my ($i, $j, $mapSum, $mapSeqLen) = ($#matrix_arrow, $#{$matrix_arrow[$#matrix_arrow]}, 0, 0);
    until ($i == 0 and $j == 0) {
        if ($matrix_arrow[$i][$j] == -1) {
            $j--;
        }
        elsif ($matrix_arrow[$i][$j] == 1) {
            $i--;
        }
        else {
            $mapSum++ if $seq1->[$i-1] eq $seq2->[$j-1];
            $i--;
            $j--;
        }
        $mapSeqLen++;
    }
	undef @matrix_score;
	undef @matrix_arrow;
    my $mapRatio = sprintf "%.3f", $mapSum/$mapSeqLen;
	my $judge = $mapRatio >= $threshold ? 1 : 0;
    return $judge, $mapRatio, 'I';
}


#function24 : output a statistical file that contains conservative or unique SSRs to the current working directory.
sub generate_cons_or_uniq_SSRs_stat_file {
	my ($fa_base_name, $out_file_suffix, $SSRs_set) = @_;
	my $out_file_name = $fa_base_name.$out_file_suffix;
	open my $OUT, '>', $out_file_name or die "Error: can't generate $out_file_name : $!";
	print $OUT "number\tid\tmotif\tmotif_length\trepeat_number\tsize\tstart\tend\tleft_fs\tleft_fs_length\tright_fs\tright_fs_length\n";
	my $count = 1;
	foreach my $fs (sort {$SSRs_set->{$a}[0] cmp $SSRs_set->{$b}[0] or $SSRs_set->{$a}[5] <=> $SSRs_set->{$b}[5]} keys %$SSRs_set){
		my ($left_fs, $right_fs) = split /-/,$fs;
		my $left_fs_len  = length $left_fs;
		my $right_fs_len = length $right_fs;
		my $out_row = join "\t", ($count, @{$SSRs_set->{$fs}}, $left_fs, $left_fs_len, $right_fs, $right_fs_len);
		print $OUT "$out_row\n";
		$count++;
	}
	close $OUT;
	undef $SSRs_set;
}


#function25 : the values of hash id_seq are connected together and stored in a scalar, and this scalar are divided equally into smaller scalars, finally, put these smaller scalars in temporary files, and return the directory name created.
sub split_id_seq_values_to_diff_file {
	my ($fileName, $dirName, $cwd, $threads_num, $idLen) = @_;
	my ($sec, $min, $hour, $day, $mon, $year, undef, undef, undef) = localtime;
	$mon  += 1;
	$year += 1900;
	my $time = "$year$mon$day-$hour$min$sec";
	my $baseName = basename $fileName;
	my $tmp_dir_name3 = $baseName."_$time.SSRMMD_tmp3";
	if (-e $tmp_dir_name3) {
		chdir $tmp_dir_name3;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name3;
	}
	mkdir $tmp_dir_name3, 0755 or die "Error: can't create '$tmp_dir_name3' directory:$!";
	chdir $cwd;
	open my $in, '<', $fileName or die "Error: can't open file '$fileName' : $!";
	chdir $dirName;
	if ($threads_num == 1) {
		my $name = $baseName.'.tmp1';
		my $outFileName = File::Spec->catfile($tmp_dir_name3, $name);
		open my $out, '>', $outFileName or die "Error: can't generate file '$outFileName' : $!";
		while (<$in>) {
			next if /\A>/;
			$_ = uc;
			s/[^A-Z]//g;
			print $out "$_\n";
		}
		close $out;
	}
	else {
		my $totalLen = 0;
		$totalLen += $idLen->{$_} foreach sort keys %$idLen;
		my $eachThreadSeqLen =  int($totalLen/$threads_num) + ($totalLen % $threads_num ? 1 : 0);
		my $nowLen = 0;
		my $count = 1;
		my $name = $baseName.".tmp$count";
		my $outFileName = File::Spec->catfile($tmp_dir_name3, $name);
		open my $out, '>', $outFileName or die "Error: can't generate file '$outFileName' : $!";
		while (<$in>) {
			next if /\A>/;
			$_ = uc;
			s/[^A-Z]//g;
			RETURNHERE:
			my $seqLen = length;
			if ($eachThreadSeqLen - $nowLen > $seqLen) {
				print $out "$_\n";
				$nowLen += $seqLen;
			}
			else {
				print $out substr($_, 0, $eachThreadSeqLen - $nowLen), "\n";
				substr($_, 0, $eachThreadSeqLen - $nowLen) = "";
				$nowLen = 0;
				$count++;
				if ($count <= $threads_num) {
					$name = $baseName.".tmp$count";
					$outFileName = File::Spec->catfile($tmp_dir_name3, $name);
					open $out, '>', $outFileName or die "Error: can't generate file '$outFileName' : $!";
				}
				my $newLength = length;
				goto RETURNHERE if $eachThreadSeqLen - $nowLen < $newLength;
				if ($_) {
					print $out "$_\n";
					$nowLen = $newLength;
				}
            }
		}
		close $out;
	}
	close $in;
	undef $idLen;
	return $tmp_dir_name3;
}


#function26 : return a array reference that contains match count of each flanking sequence compared with genome sequence.
sub get_fs_match_count {
	my ($file, $new_fs_set, $fs_len) = @_;
	open my $in, '<', $file or die "Error: can't open file '$file' : $!";
	undef $file;
	my $match_count;
	@$match_count = (0) x @$new_fs_set;
	my ($short_seq, $shseq_len) = ('', 0);
	while (<$in>) {
		s/[\r\n]+//;
		$short_seq .= $_;
		$shseq_len += length;
		if ($shseq_len > 10_000_000) {
			sldingWindows($new_fs_set, $fs_len, \$short_seq, $shseq_len, $match_count);
			($short_seq, $shseq_len) = ('', 0);
		}
	}
	close $in;
	if ($shseq_len ) {
		sldingWindows($new_fs_set, $fs_len, \$short_seq, $shseq_len, $match_count);
	}
	undef $new_fs_set;
	undef $short_seq;
	return $match_count;
}


#function27 : check uniqueness of flanking sequences by slding window.
sub sldingWindows {
	my ($new_fs_set, $fs_len, $short_seq, $shseq_len, $match_count) = @_;
	foreach my $first_start_loci (0 .. ($fs_len-1)){ #the number of hash %seq_count generated is the same as the length of the flanking sequence
		my %seq_count = ();
		my $start_loci = $first_start_loci;
		foreach my $num (1 .. int( ($shseq_len - $first_start_loci) / $fs_len )){ #the number of cuts when the genome sequence is cut to the same length as the flanking sequence
			my $seq = substr $$short_seq, $start_loci, $fs_len;
			$seq_count{$seq}++;
			$start_loci += $fs_len;
		}
		my $index = -1;
		foreach my $i (0 .. $#{$new_fs_set}){ #each flanking sequence is aligned with the genome sequence
			$index++;
			my $count = 0;
			foreach my $each_fs (@{$new_fs_set->[$i]}){ #flanking sequence is divided into left and right
				next unless $seq_count{$each_fs};
				if($seq_count{$each_fs} > 1){ #flanking sequence on either side is matched to more than one position
					$count = 3;
					last;
				}
				$count += 1;
			}
			$match_count->[$index] += $count;
		}
	}
	undef $new_fs_set;
	undef $short_seq;
	#notice : By doing this, the flank sequence(140bp) that is not enough to 
	####### : set the value(150bp) will default to mismatch. Obviously, this 
	####### : kind of situation rarely happen. However, we have taken remedial 
	####### : measure: if the final match count is less than or equal to 2, 
	####### : the flanking sequence will be regarded as unique.
}


#function28 : use multithreads and call function get_fs_match_count to keep the unique SSRs(the flanking sequence is unique) in the genome, and return a hash reference.
sub multithreads_keep_unique_SSRs {
	my ($tmp_dir_name3, $fs_con_SSRs_file, $fs_len) = @_;
	die "Error: no find file '$fs_con_SSRs_file' : $!" unless -e $fs_con_SSRs_file;
	my $fs_con_SSRs = retrieve $fs_con_SSRs_file;
	my $fs_set = [ sort {$fs_con_SSRs->{$a}[0] cmp $fs_con_SSRs->{$b}[0] or $fs_con_SSRs->{$a}[5] <=> $fs_con_SSRs->{$b}[5]} keys %$fs_con_SSRs ];
	undef $fs_con_SSRs;
	my $new_fs_set;
	foreach my $fs (@$fs_set){
		my @fs = split /-/, $fs;
		push @{$new_fs_set}, [@fs];
	}
	undef $fs_set;
	unless(defined $new_fs_set){
		if(-e $tmp_dir_name3){
			chdir $tmp_dir_name3;
			unlink glob '*';
			chdir '..';
			rmdir $tmp_dir_name3;
		}
		undef $tmp_dir_name3;
		unlink $fs_con_SSRs_file;
		undef $fs_con_SSRs_file;
		goto LOOP;
	}
	my @all_count = ();
	my @files = glob "$tmp_dir_name3/*";
	if (@files == 1) {
		my $match_count = get_fs_match_count($files[0], $new_fs_set, $fs_len);
		undef $new_fs_set;
		push @all_count, $match_count;
	}
	elsif (@files >= 2) {
		foreach my $file (@files){
			my $thr = threads->create(\&get_fs_match_count, $file, $new_fs_set, $fs_len);
		}
		undef $new_fs_set;
		while(threads->list()){
    		foreach my $thr (threads->list(threads::joinable)){
        		my $match_count = $thr->join();
				push @all_count, $match_count;
    		}
		}
	}
	else {
		die "Error: can't find any file in tmp directory '$tmp_dir_name3' : $!";
	}
	undef $new_fs_set;
	if (-e $tmp_dir_name3) {
		chdir $tmp_dir_name3;
		unlink glob '*';
		chdir '..';
		rmdir $tmp_dir_name3;
	}
	undef $tmp_dir_name3;
	my @match_count = ();
	foreach my $i (0 .. $#{$all_count[0]}){
		my $count = 0;
		$count += $all_count[$_][$i] foreach 0 .. $#all_count;
		push @match_count, $count;
	}
	undef @all_count;
	my %fs_uni_SSRs;
	my $index = 0;
	$fs_con_SSRs = RetrieveFromFile($fs_con_SSRs_file);
	undef $fs_con_SSRs_file;
	foreach my $fs (sort {$fs_con_SSRs->{$a}[0] cmp $fs_con_SSRs->{$b}[0] or $fs_con_SSRs->{$a}[5] <=> $fs_con_SSRs->{$b}[5]} keys %$fs_con_SSRs){
		$fs_uni_SSRs{$fs} = $fs_con_SSRs->{$fs} if $match_count[$index] <= 2; #SSR with flanking sequence matching count of 0, 1, 2 is retained
		$index++;
	}
	undef @match_count;
	undef $fs_con_SSRs;
	LOOP:
	return \%fs_uni_SSRs;
}


#function29 : output a statistical table that contains polymorphic SSRs candidate to the current working directory.
sub final_comparison_table {
	my ($fs_uni_SSRs1, $fs_uni_SSRs2, $con_compared, $fa_base_name1, $fa_base_name2, $conserMap) = @_;
	my $out_file_name = $fa_base_name1.'-and-'.$fa_base_name2.'.compare';
	open my $OUT,'>',$out_file_name or die "Error: can't generate $out_file_name : $!";
	print $OUT "number\tfasta1_id\tfasta1_motif\tfasta1_repeat_number\tfasta1_start\tfasta1_end\tfasta2_id\tfasta2_motif\tfasta2_repeat_number\tfasta2_start\tfasta2_end\tfasta1_left_fs\tfasta1_left_fs_length\tfasta2_left_distance(LD)\tfasta2_left_identity(NW)\tfasta1_right_fs\tfasta1_right_fs_length\tfasta2_right_distance(LD)\tfasta2_right_identity(NW)\tpolymorphism\n";
	my $count = 1;
	foreach my $fs (sort {$con_compared->{$a}[0] cmp $con_compared->{$b}[0] or $con_compared->{$a}[5] <=> $con_compared->{$b}[5]} keys %$con_compared){
		my $fs2 = $con_compared->{$fs}[-5] eq 'same' ? $fs : $con_compared->{$fs}[-5];
		if(exists $fs_uni_SSRs1->{$fs} and exists $fs_uni_SSRs2->{$fs2}){
			my ($Lresult, $Lwhich, $Rresult, $Rwhich) = @{$con_compared->{$fs}}[-4 .. -1];
			my $Ldistance = $Lwhich eq 'D' ? $Lresult : 'NA';
			my $Lidentity = $Lwhich eq 'I' ? $Lresult : 'NA';
			my $Rdistance = $Rwhich eq 'D' ? $Rresult : 'NA';
			my $Ridentity = $Rwhich eq 'I' ? $Rresult : 'NA';
			if ($Lresult eq 'NA' and $Lwhich eq 'NA' and $Rresult eq 'NA' and $Rwhich eq 'NA') {
				if ($conserMap eq 'NW') {
					$Lidentity = '1.000';
					$Ridentity = '1.000';
				}
				else {
					$Ldistance = '0.000';
					$Rdistance = '0.000';
				}
			}
			my ($left_fs, $right_fs) = split /-/, $fs;
			my $left_fs_len  = length $left_fs;
			my $right_fs_len = length $right_fs;
			my $judge = $con_compared->{$fs}[3] != $con_compared->{$fs}[10] ? 'yes' : 'no';
			my $out_row = join "\t", ($count, @{$con_compared->{$fs}}[0,1,3,5,6, 7,8,10,12,13], $left_fs, $left_fs_len, $Ldistance, $Lidentity, $right_fs, $right_fs_len, $Rdistance, $Ridentity, $judge);
			print $OUT "$out_row\n";
			$count++;
		}
	}
	close $OUT;
	undef $fs_uni_SSRs1;
	undef $fs_uni_SSRs2;
	undef $con_compared;
}
