#! /usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;


##get options from command line into variables and set default values
my ($thread_num, $buffer_size, $qualSys, $Help, $P, $q);
GetOptions(
	"t:i"=>\$thread_num,
	"m:i"=>\$buffer_size,
	"Q:i"=>\$qualSys,
	"p:s"=>\$P,
	"q:s"=>\$q,
	"h"=>\$Help
);
die "Version: 5.1
Usage:$0 <lane.lst> <lib.lst> [maxjob, default 10]
        -t <int>    thread number for each filter_data job,default 8
	-m <int>    the reads pair number in buffer,default 2000000
	-Q <int>    quality system: 64 or 33 (default: 33)
	-p <str>    set the qsub -P,default=rdtest
        -q <str>    set the qsub -q,default=bc.q
	-h          output help information to screen\n" if(@ARGV<2||$Help);
my $lane_lst = shift;
my $lib_lst = shift;
my $Maxjob = shift;

$Maxjob ||= 10;
$thread_num ||= 8;
$buffer_size ||= 2000000;
$qualSys ||= 33;
$P ||= "rdtest";
$q ||= "bc.q";

my $file_name = `basename $lane_lst`;
my $dir_name=`pwd`;
chomp $file_name;
chomp $dir_name;

my %Lib;
my $vf="";
open IN,$lib_lst or die "$!";
while(<IN>){
	if (/^(\S+)\s+(\d+)/){
		$Lib{$1}=$2;
		mkdir($1) unless (-d $1);
	}else{
		next;
	}
}
close IN;

open OUT1, ">$file_name.filter.sh" or die "$!";
open OUT2, ">$file_name.dup.sh" or die "$!";
open OUT3, ">$file_name.stat.sh" or die "$!";
open OUT4, ">$file_name.checkfq.sh" or die "$!";

if ( -e "$file_name.stat.xls" ){
	my $tag = time;
	print STDERR "$file_name.stat.xls  already exists! \n";
	print STDERR "mv  $file_name.stat.xls $file_name.stat.xls.$tag \n";
	`mv  $file_name.stat.xls $file_name.stat.xls.$tag`;
}

open IN,$lane_lst or die "$!";

my %raw_path;
my $stat_output = '';
my $num = 0;
while(<IN>){
	next if !/_1\.fq/;
	chomp;
	$num++;
	
	my ($f1,$start1,$end1,$B_cutoff) = split(/\s+/);
	$start1 ||= 0;
	$end1 ||= 0;
	$B_cutoff ||= 40;
	my $raw_path0=`dirname $f1`;
	chomp $raw_path0;
	my $raw_path1=`dirname $raw_path0`;
	chomp $raw_path1;
	my $raw_path2=`dirname $raw_path1`;
	chomp $raw_path2;
	$raw_path{"$raw_path2"}=1;
	
	my $line2 = <IN>;
	chomp $line2;
	my ($f2,$start2,$end2,$N_num) = split(/\s+/,$line2);
	$start2 ||= 0;
	$end2 ||= 0;
	#$N_num ||= 10;
	if(not defined $N_num or $N_num eq "")
	{$N_num ||= 10;}
	
	my $name = `basename $f1`;
	chomp $name;
	my $name2 = `basename $f2`;
	chomp $name2;
    my $prefix=$1 if ($name=~/(.*)_1/);
	my $lib;
	if ( $name =~/L\d+_([^_]+)_1\.fq/ ){
		$lib = $1;
	}else{
		die;
	}
	next if not exists $Lib{$lib};

	print OUT1 "$Bin/filter_data_parallel -t $thread_num -m $buffer_size -q $qualSys -y -z -w $N_num -B $B_cutoff -l $Lib{$lib} -a $start1 -b $end1 -c $start2 -d $end2 $f1 $f2 ./$lib/$name.reads.stat ./$lib/$name.clean.tmp  ./$lib/$name2.clean.tmp && mv ./$lib/$name.clean.tmp ./$lib/$name.clean && mv ./$lib/$name2.clean.tmp ./$lib/$name2.clean && echo OK \n" ;
	
	print OUT2 "$Bin/duplication ./$lib/$name.clean ./$lib/$name2.clean ./$lib/$name.clean.dup.clean.gz ./$lib/$name2.clean.dup.clean.gz ./$lib/$name.clean.dup.stat && rm ./$lib/$name.clean ./$lib/$name2.clean && echo OK \n";
	#print OUT2 "$Bin/duplication ./$lib/$name.clean ./$lib/$name2.clean ./$lib/$name.clean.dup.clean ./$lib/$name2.clean.dup.clean ./$lib/$name.clean.dup.stat && rm ./$lib/$name.clean ./$lib/$name2.clean && echo OK \n";
#	my $parameter="$start1"."_"."$end1"."_"."$start2"."_"."$end2"."_"."$B_cutoff"."_"."$N_num";
#	print OUT3 "$Bin/stat.pl $lib_lst ./$lib/$name.reads.stat ./$lib/$name.clean.dup.stat $parameter >>$file_name.stat.xls \n";

        if($qualSys == 33){
	  print OUT4 "$Bin/checkfastq_s33 $dir_name/$lib/$name.clean.dup.clean.gz > $dir_name/$lib/$name.clean.dup.clean.gz.check\n$Bin/checkfastq_s33 $dir_name/$lib/$name2.clean.dup.clean.gz > $dir_name/$lib/$name2.clean.dup.clean.gz.check\nperl $Bin/distribute_fqcheck.pl $dir_name/$lib/$name.clean.dup.clean.gz.check $dir_name/$lib/$name2.clean.dup.clean.gz.check -o $dir_name/$lib/$prefix && echo OK\n";
	}elsif($qualSys == 64){
	  print OUT4 "$Bin/checkfastq_s64 $dir_name/$lib/$name.clean.dup.clean.gz > $dir_name/$lib/$name.clean.dup.clean.gz.check\n$Bin/checkfastq_s64 $dir_name/$lib/$name2.clean.dup.clean.gz > $dir_name/$lib/$name2.clean.dup.clean.gz.check\nperl $Bin/distribute_fqcheck.pl $dir_name/$lib/$name.clean.dup.clean.gz.check $dir_name/$lib/$name2.clean.dup.clean.gz.check -o $dir_name/$lib/$prefix && echo OK\n";
	}

	my $parameter="$start1"."_"."$end1"."_"."$start2"."_"."$end2"."_"."$B_cutoff"."_"."$N_num";
	print OUT3 "$Bin/stat.pl $lib_lst ./$lib/$name.reads.stat ./$lib/$name.clean.dup.stat $parameter >>$file_name.stat.xls \n";


	if($Lib{$lib}<1000){$vf="7G";}
	else{$vf="5G";} 

}

close IN;
close OUT1;
close OUT2;
close OUT3;
close OUT4;

print STDERR "Runing low quality filtering ... \n";
`perl $Bin/qsub-sge.pl -maxjob $Maxjob -resource="vf=3.5G -q $q -P $P" -reqsub $file_name.filter.sh`;
print STDERR "Finished. \n";
print STDERR "Runing duplicate filtering ... \n";
`perl $Bin/qsub-sge.pl -maxjob $Maxjob -resource="vf=$vf -q $q -P $P" -reqsub $file_name.dup.sh`;
print STDERR "Finished. \n";
print STDERR "Checking fastq reads at back ground ...\n";
`perl $Bin/qsub-sge.pl -maxjob $Maxjob -resource="vf=2G -q $q -P $P" -reqsub  -lines 3  $file_name.checkfq.sh`;
print STDERR "Finished. \n";
print STDERR "Runing stat. ... \n";
`sh $file_name.stat.sh`;
print STDERR "Finished. \n";

open OUT4,">$file_name.stat.new.sh";
foreach my $aa(sort keys %raw_path){
   print OUT4 "perl $Bin/stat_process.pl $file_name.stat.xls $aa $dir_name >>$file_name.stat.tmp.xls\n";
}
print OUT4 "cat $Bin/stat_title.txt $file_name.stat.tmp.xls >$file_name.stat.new.xls; rm $file_name.stat.tmp.xls\n";
close OUT4;

print STDERR "Runing stat new ... \n";
`sh $file_name.stat.new.sh`;
print STDERR "All finished. \n";

