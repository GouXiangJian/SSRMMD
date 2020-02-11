#!/usr/bin/env perl

#date    : 2020-02-10
#author1 : Xiang-jian Gou (xjgou@stu.sicau.edu.cn)
#author2 : Hao-ran Shi (542561234@qq.com)

#You can also set more primer design conditions, in script line 249-262.

#load modules
use strict;
use warnings;
use Cwd qw/abs_path/;
use Getopt::Long;

#clear buffer
$| = 1;

#record version information
my $VERSION = 'connectorToPrimer3 v1.0';

#get platform information
my $platform = $^O; #MSWin32 linux

#get path of connectorToPrimer3
my $path = abs_path $0;
$path =~ s/\\/\//g;
$path =~ s/\/[\w\.]+\z//;

#set default options
my $inFile = '';
my $outFile = 'OUTPUT';
my $inFileSource = 1;
my $mode = 'USE'; #IN OUT USE
my $filter = 0;
my $returnPrimerSum = 1;
my $GC = '20,80';
my $TM = '57,60,63';
my $TMdiff = 2;
my $primerSize = '18,20,27';
my $produceSize = '100-200';
my $version;
my $help;

#get options from command line
GetOptions(
        'input=s'           =>  \$inFile,
        'output=s'          =>  \$outFile,
        'source=i'          =>  \$inFileSource,
        'mode=s'            =>  \$mode,
        'filter=i'          =>  \$filter,
        'path=s'            =>  \$path,
        'return=i'          =>  \$returnPrimerSum,
        'gc=s'              =>  \$GC,
        'tm=s'              =>  \$TM,
        'diff=i'            =>  \$TMdiff,
        'primerSize|es=s'   =>  \$primerSize,
        'produceSize|us=s'  =>  \$produceSize,
        'version+'          =>  \$version,
        'help+'             =>  \$help,
);

#describe program information
my $usage = <<__GUIDE__;
###############################################################################
Name:
  connectorToPrimer3 - Provide the connection between SSRMMD and Primer3

Author:
  Xiang-jian Gou (xjgou\@stu.sicau.edu.cn)
  Hao-ran Shi (542561234\@qq.com)

Usage:
  perl connectorToPrimer3.pl opt1 <value1> opt2 <value2> ... optN <valueN>

Options:
  -i  | -input  <STR> : provide an input file name, just accept output file
                        from SSRMMD or Primer3. (must be provided)

  -o  | -output <STR> : create a file for storing output. (default: OUTPUT)

  -s  | -source <INT> : the source of input(option '-i') file. (default: 1)
                           1 : from SSRMMD when just mining SSR loci 
                           2 : from SSRMMD when mining polymorphic loci
                           3 : from Primer3

  -m  | -mode   <STR> : select a mode to run connectorToPrimer3. (default: USE)
                           IN  : generate Primer3(P3) input
                           OUT : parse Primer3 output
                           USE : generate P3 input, run P3, and parse P3 output

  -f  | -filter <INT> : if option -s = 2, filter input file. (default: 0)
                           0 : no filter
                           1 : just deal polymorphic loci (col 20 = 'yes')
                           2 : just deal distance(LD) = 0 (col 14,18 = '0.000')
                           3 : just deal identity(NW) = 1 (col 15,19 = '1.000')
                           4 : meet both 1 and 2
                           5 : meet both 1 and 3

  -pa | -path   <STR> : if option -m = USE, specify path of primer3_core-linux
                        or primer3_core-win.exe. (default: the path same as
                        connectorToPrimer3.pl)

<Here are some options of primer design>:

  -r  | -return <INT> : number of primers each SSR want to design. (default: 1)

  -g  | -gc     <STR> : set minimum and maximum GC content. (default: 20,80)

  -t  | -tm     <STR> : set min, opt and max TM value. (default: 57,60,63)

  -d  | -diff   <INT> : max difference in TM between primer pair. (default: 2)

  -es | -primerSize  <STR> : set min,opt,max primer size. (default: 18,20,27)

  -us | -produceSize <STR> : set produce size range. (default: 100-200)

  -v  | -version           : show the version information.

  -h  | -help              : show the help information.
###############################################################################

__GUIDE__

#parse some options
my ($minGC, $maxGC) = $GC =~ /\A(\d+),(\d+)\z/;
my ($minTM, $optTM, $maxTM) = $TM =~ /\A(\d+),(\d+),(\d+)\z/;
my ($minPriSize, $optPriSize, $maxPriSize) = $primerSize =~ /\A(\d+),(\d+),(\d+)\z/;

#check the options
die "$VERSION\n" if $version;
die $usage if $help;
die "Error: option '-i'  must be provided !\n" if ! $inFile;
die "Error: option '-s'  just be set to 1, 2 or 3 !\n" if $inFileSource != 1 and $inFileSource != 2 and $inFileSource != 3;
die "Error: option '-m'  just be set to 'IN', 'OUT' or 'USE' !\n" if $mode ne 'IN' and $mode ne 'OUT' and $mode ne 'USE';
die "Error: option '-f'  just be set to 0, 1, 2, 3, 4 or 5 !\n" if $filter != 0 and $filter != 1 and $filter != 2 and $filter != 3 and $filter != 4 and $filter != 5;
die "Error: option '-r'  must be an integer bigger than 0 !\n" if $returnPrimerSum <= 0 or $returnPrimerSum =~ /\.\d*[1-9]+/;
die "Error: option '-g'  have a problem !\n" if (! defined($minGC)) or (! defined($maxGC));
die "Error: option '-t'  have a problem !\n" if (! defined($minTM)) or (! defined($optTM)) or (! defined($maxTM));
die "Error: option '-d'  must be an integer bigger than 0 !\n" if $TMdiff < 0 or $TMdiff =~ /\.\d*[1-9]+/;
die "Error: option '-es' have a problem !\n" if (! defined($minPriSize)) or (! defined($optPriSize)) or (! defined($maxPriSize));
die "Error: option '-us' have a problem !\n" if $produceSize !~ /\A(\d+)-(\d+)\z/;
die "Error: when option '-m' = 'IN',  option '-s' can't be set to 3 !\n" if $mode eq 'IN' and $inFileSource == 3;
die "Error: when option '-m' = 'OUT', option '-s' can't be set to 1 or 2 !\n" if ($mode eq 'OUT' and $inFileSource == 1) or ($mode eq 'OUT' and $inFileSource == 2);
die "Error: when option '-m' = 'USE', option '-s' can't be set to 3 !\n" if $mode eq 'USE' and $inFileSource == 3;



#######################
#start main program ...
#######################



if ($mode eq 'IN') {
    parseSSRMMD($outFile);
}
elsif ($mode eq 'OUT') {
    parsePrimer3($inFile);
}
else {
    my ($sec, $min, $hour, $day, $mon, $year, undef, undef, undef) = localtime;
	$mon  += 1;
	$year += 1900;
	my $time = "$year$mon$day-$hour$min$sec";
	my $fileName1 = "connectorToPrimer3.$time.tmp1";
    my $fileName2 = "connectorToPrimer3.$time.tmp2";
    parseSSRMMD($fileName1);
    chop $path if $path =~ /\/\z/;
    if ($platform eq 'MSWin32') {
        #primer3 in windows is v2.5.0
        #return 0 imply that 'system' run right !
        system "$path/primer3_core-win.exe $fileName1 --output=$fileName2" and die "Error: can't run $path/primer3_core-win.exe : $!";
    }
    elsif ($platform eq 'linux') {
        #primer3 in linux is v2.5.0
        system "chmod 755 $path/primer3_core-linux" and die "Error: failed to modify file $path/primer3_core-linux permissions : $!";
        system "$path/primer3_core-linux $fileName1 --output=$fileName2" and die "Error: can't run $path/primer3_core-linux : $!";
    }
    else {
        warn "Warning: your platform isn't windows or linux, try to run primer3 by primer3_core-linux !\n";
        system "chmod 755 $path/primer3_core-linux" and die "Error: failed to modify file $path/primer3_core-linux permissions : $!";
        system "$path/primer3_core-linux $fileName1 --output=$fileName2" and die "Error: can't run $path/primer3_core-linux : $!";
    }
    unlink $fileName1 if -e $fileName1;
    parsePrimer3($fileName2);
    unlink $fileName2 if -e $fileName2;
}



#####################
#end main program ...
#####################



#parse output of SSRMMD
sub parseSSRMMD {
    my $outName = shift;
    open my $in, '<', $inFile or die "Error: can't open file '$inFile' : $!";
    my $first = <$in>;
    die "Error: the input file '$inFile' have a wrong format !\n" if (! defined $first) or $first !~ /\Anumber/;
    my $sum = split /\t/, $first;
    die "Error: when option '-s' = 1, the input file '$inFile' have a wrong format !\n" if $sum != 12 and $inFileSource == 1;
    die "Error: when option '-s' = 2, the input file '$inFile' have a wrong format !\n" if $sum != 20 and $inFileSource == 2;
    open my $out, '>', $outName or die "Error: can't write file '$outName' : $!";
    while (<$in>) {
        s/[\r\n]+//;
        my @row = split /\t/;
        if ($inFileSource == 2 and $filter != 0) {
            my $polyJudge = $row[-1] ne 'yes';
            my $ldJudge = $row[13] ne '0.000' || $row[17] ne '0.000';
            my $nwJudge = $row[14] ne '1.000' || $row[18] ne '1.000';
            if ($filter == 1) {
                next if $polyJudge;
            }
            elsif ($filter == 2) {
                next if $ldJudge;
            }
            elsif ($filter == 3) {
                next if $nwJudge;
            }
            elsif ($filter == 4) {
                next if $polyJudge or $ldJudge;
            }
            elsif ($filter == 5) {
                next if $polyJudge or $nwJudge;
            }
            else {
                die "Error : option '-f' just be set to 0, 1, 2, 3, 4 or 5 !\n";
            }
        }
        my $seqId = $row[0];
        my ($seq, $start, $length);
        if ($inFileSource == 1) {
            $seq = $row[8].($row[2]x$row[4]).$row[10];
            $start = $row[9]+1;
            $length = $row[5];
        }
        elsif ($inFileSource == 2) {
            $seq = $row[11].($row[2]x$row[3]).$row[15];
            $start = $row[12]+1;
            $length = length($row[2])*$row[3];
        }
        else {
            die "Error: program run have wrong (option '-s' = $inFileSource) !\n";
        }
        print $out <<__P3IN__
SEQUENCE_ID=$seqId
SEQUENCE_TEMPLATE=$seq
SEQUENCE_TARGET=$start,$length
PRIMER_NUM_RETURN=$returnPrimerSum
PRIMER_MIN_GC=$minGC
PRIMER_MAX_GC=$maxGC
PRIMER_MIN_TM=$minTM
PRIMER_MAX_TM=$maxTM
PRIMER_OPT_TM=$optTM
PRIMER_PAIR_MAX_DIFF_TM=$TMdiff
PRIMER_MIN_SIZE=$minPriSize
PRIMER_MAX_SIZE=$maxPriSize
PRIMER_OPT_SIZE=$optPriSize
PRIMER_PRODUCT_SIZE_RANGE=$produceSize
=
__P3IN__
    }
    close $in;
    close $out;
}

#parse output of Primer3
sub parsePrimer3 {
    my $inName = shift;
    open my $in, '<', $inName or die "Error: can't open file '$inName' : $!";
    my $first = <$in>;
    die "Error: the input file '$inName' have a wrong format !\n" if (! defined $first) or $first !~ /\ASEQUENCE_ID=\d+/;
    seek $in, 0, 0;
    local $/ = "=\n";
    open my $out, '>', $outFile or die "Error: can't write file '$outFile' : $!";
    print $out join("\t", qw/id forward_seq forward_length forward_TM forward_GC reverse_seq reverse_length reverse_TM reverse_GC produce_size/), "\n";
    while (<$in>) {
        my ($total) = /PRIMER_PAIR_NUM_RETURNED=(\d+)/;
        next if (! defined $total) or $total == 0;
        my ($id) = /SEQUENCE_ID=(\w+)/;
        foreach my $i (0 .. $total-1) {
            my ($leftSeq)   = /PRIMER_LEFT_${i}_SEQUENCE=(\w+)/;
            my $leftSeqLen  = length $leftSeq;
            my ($rightSeq)  = /PRIMER_RIGHT_${i}_SEQUENCE=(\w+)/;
            my $rightSeqLen = length $rightSeq;
            my ($leftTM)  = /PRIMER_LEFT_${i}_TM=([\d\.]+)/;
            my ($rightTM) = /PRIMER_RIGHT_${i}_TM=([\d\.]+)/;
            my ($leftGC)  = /PRIMER_LEFT_${i}_GC_PERCENT=([\d\.]+)/;
            my ($rightGC) = /PRIMER_RIGHT_${i}_GC_PERCENT=([\d\.]+)/;
            my ($realProduceSize) = /PRIMER_PAIR_${i}_PRODUCT_SIZE=(\d+)/;
            my $realId = "$id.".($i+1);
            print $out join("\t", $realId, $leftSeq, $leftSeqLen, $leftTM, $leftGC, $rightSeq, $rightSeqLen, $rightTM, $rightGC, $realProduceSize), "\n";
        }
    }
    close $in;
    close $out;
}
