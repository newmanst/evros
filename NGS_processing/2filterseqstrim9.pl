#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw( gettimeofday);
use Cwd;
use Exporter;
use Carp;
use List::Util ();
use File::Copy;
use threads;
#    arguments:
#    perl filterseqs.pl          52			 52                       100                    20	     	  1
#    ____    this           min  length     max length        minimum % allowed   minimum Q allowed     check quality

my $help="something has gone wrong. help";
my $OS_ERROR="OS error!";

#if help is requested, output help string, and end program
if(exists($ARGV[0])&&$ARGV[0]eq"help"){	print $help;	die "goodbye";}

######## editable params #############
##### edit based on application ######
	my $lmin=40;
	my $lmax=65;
######################################
######################################



######### check inputs ###############
#default params
my $minlength=52;
my $maxlength=52;
my $percentallowed=99;
my $minimumQ=20;
my $checkquality=0;


#replace default params with defined
if(exists($ARGV[0])){$minlength=$ARGV[0];}
if(exists($ARGV[1])){$maxlength=$ARGV[1];}
if(exists($ARGV[2])){$percentallowed=$ARGV[2];}
if(exists($ARGV[3])){$minimumQ=$ARGV[3];}
if(exists($ARGV[4])){$checkquality=$ARGV[4];}

######################################


############## read current directory #####################
my $some_dir=getcwd;
my @directories=();
my @foldernames=();
opendir(my $dh, $some_dir) || die "Can't open $some_dir: $!";
my $count=0;
####line 11###
while (readdir $dh) {
	if(index($_,".")==-1){
		$directories[$count]=$some_dir."/".$_;
		$foldernames[$count]=$_;
		print "".($count+1).") ".$_."\n";
		$count++;
	}
}
print "\nWhich directory would you like to choose?";
my $choice=<STDIN>;
print "\n";
closedir $dh;
###########################################################

################ output list of files to test #############
my @filenames=getFiles($some_dir,$foldernames[$choice-1]);
my $in_location="".$foldernames[$choice-1]."/Extracted/";
my $out_location="".$foldernames[$choice-1]."/Filtered/";

#open(my $out_file,'>',$out_location."FastqProcesser2_output.txt");

my $new_directory=$foldernames[$choice-1]."/Filtered/";
unless(-e $new_directory or mkdir($new_directory, 0755)) {
        die "Unable to create $new_directory\n";
}

############## start analysis #######################
#create file to store run parameters
open(my $params, '>',$foldernames[$choice-1]."/FilterOutput.txt")||die "could not open ".$foldernames[$choice-1]."/FilterOutput.txt";
print $params "Min Length = $minlength\nMax Length = $maxlength\nAllowing sequences with at least $percentallowed percent of bases with Q>=$minimumQ\n";
#print "Min Length = $minlength\nMax Length = $maxlength\nAllowing sequences with at least $percentallowed percent of bases with Q>=$minimumQ\n";


############## start analysis #######################
my @thr=();
for(my $q=0;$q<scalar(@filenames);$q++)#31;$q<=31;$q++)#
{


	my $file = $filenames[$q];
	open(INPUT,$in_location.$file) || die("\nERROR: Invalid File ".$in_location.$file); 
	print $in_location.$file."\n";
	my $line_pos=0;
	my $id1="";
	my $seq="";
	my $id2="";
	my $Q="";
	my $count=0;
	my $qcount=0;
	my %unique=();
	my @splitfile=split("_",$file);

	#rename output file --> 'XX##_R1.fastq'
	my $newfilename=$splitfile[0]."_".$splitfile[3].".fastq";

	######################### REMOVE IF COPIED ################################
	######## ************************************************ #################
	######## ************************************************ #################
	######## only to be used if there is no zero in files 1-9 #################
	######## ************************************************ #################
	######## ************************************************ #################
    #if(length($newfilename)<13){
		#print "did you remove the name check from the previous iteration of the program??\n";
    #	$newfilename="BW0".substr($newfilename,2);
    #}
	######## ************************************************ #################
	######## ************************************************ #################
	######## ************************************************ #################
	######## ************************************************ #################

	#count lengths=48,4950,51,52,53,54,55
	my @lengthcount=();
	my @lengthcount_unique=();
	my @chimera_count=();
    if($q == 0){print $params "Filename\tTotalReads\tPassingQ";}
	for(my $i=$lmin;$i<=$lmax;$i++){
		$lengthcount[$i-$lmin]=0;
        if($q == 0){print $params "\tL=".$i;}
    }
    if($q == 0){print $params "\n";}
    my $scount=0;
	open(my $nf,'>',$out_location.$newfilename)|| die("\nERROR: Invalid File ".$out_location.$newfilename);
	while(my $in_line = <INPUT>) 
	{	if($scount<10000000){
		$in_line=~s/\n//g;
		if($line_pos==0){
			$id1=$in_line;
		}
		if($line_pos==1){
			$seq=$in_line;
		}
		if($line_pos==2){
			$id2=$in_line;
		}
		$line_pos++;
		if($line_pos==4){
		
			$Q=$in_line;
		##################################################
		##################################################
		###### for removing the last 9 bases #############
		##################################################
            if(length($seq)>=52){
                $seq=substr($seq,0,52);
                    
                $Q=substr($Q,0,52);
            }
		##################################################
			$line_pos=0;
			my $l=length($seq);
			if($l>=$lmin&&$l<=$lmax){
				$lengthcount[$l-$lmin]+=1;
			}
			if($checkquality)
			{
					$scount++;
				if(goodquality($percentallowed,$minimumQ,$minlength,$maxlength,$Q)==1){
					$qcount++;
					print $nf $id1."\n".$seq."\n".$id2."\n".$Q."\n";
				}
			}
			else{
				$qcount++;
				print $nf $id1."\n".$seq."\n".$id2."\n".$Q."\n";
			}
			$count++;
		}
		}
	}
	
	close($nf);
	close(INPUT);
    print $params "".getfilename($newfilename)."\t".$count."\t".$qcount;
    for(my $i=$lmin;$i<=$lmax;$i++){
        print $params "\t".$lengthcount[$i-$lmin];
    }
    print $params "\n";
	print "Success. Converted ".getfilename($file)."-->".getfilename($newfilename)."\t".$count."\t".$qcount."\n";
}

close($params);


#args: counter,start directory
sub getDirectories{
	my $shifter=shift;
	my $dir=shift;
    opendir(my $dh, $dir) || die "shit";
    while (readdir $dh) {
	if(index($_,".")==-1){
	    for(my $i=0;$i<$shifter;$i++){print"\t";}
            print "$dir/$_\n";
	    getDirectories($shifter+1,"$dir/$_");
	}
    }
    closedir $dh;
}

#args: counter,start subfolder
sub getFiles{

	my $dir = shift;
	my $folder = shift;
	my @files=();
	my $count=0;
	opendir(my $dh, $dir."/".$folder."/Extracted") || die "shit";
	while (readdir $dh) {
		if($_ ne "." && $_ ne ".."){
			print"\t";
        		print "$_\n";
        		$files[$count]=$_;
        		$count++;
		}
    	}
	closedir $dh;
	return @files;
}
#counter, origin, destination, start directory
sub copyDirectories{ #level 2-->0
	my $shifter=shift;
	my $origin=shift;
	my $dest=shift;
	my $dir=shift;
	my $flag=0;
	if($shifter==1){
		print $dir."\n";
	}
	opendir(my $dh, $dir) || die "shit";
	while (readdir $dh) {
		for(my $i=0;$i<$shifter;$i++){
			print"\t";
		}
        	print "$dir/$_\n";
		if(index($_,".")==-1){
			copyDirectories($shifter+1,$origin,$dest,"$dir/$_");
		}else{
			#print $shifter."\t".$origin."\n";
			if($shifter==$origin && $_ ne "." && $_ ne ".."){
				my $sourcefile=truncDirectory(5,$dir."/".$_);
				my $destdir=truncDirectory(5,moveDirectory($origin-$dest,$dir."/".$_)."/Extracted/")."/";
				system("7zG x ".$sourcefile." -o".$destdir);#."\n";
				for(my $i=0;$i<$shifter;$i++){
					print"\t";
				}
				print "    --> ";
				print " successfully extracted\n";# $sourcefile to $destdir\n";
			}
		}
	}
	closedir $dh;
}

#args: number of moves upwards in directory, subdirectory
sub moveDirectory{
	my $delta=shift;
	my $dir=shift;
	my @temp=split("/",$dir);
	for(my $i=1;$i<=$delta;$i++){
		$temp[scalar(@temp)-$i]="";
	}
	my $output=join("/",@temp);
	return substr($output,0,length($output)-$delta);
}

#args: number of prior truncations, full directory
sub truncDirectory{
	my $delta=shift;
	my $dir=shift;
	my @temp=split("/",$dir);
	for(my $i=0;$i<$delta;$i++){
		$temp[$i]="";
	}
	my $output=join("/",@temp);
	return substr($output,$delta,length($output)-$delta);

}

#args: name of subfolder (from choice above)
sub unzipAll{
	my $subfolder=shift;
	system("7za a -$subfolder *.fastq");
}
sub goodquality{
    
    my $percentallowed=$_[0]; #98
    my $Qfloor=$_[1]; #18
    my $minLength=$_[2]; #51
    my $maxLength=$_[3]; #53
    splice @_, 0, 4;

    my $totalitems=0;
    my $passingitems=0;
    
    my @str=split("",$_[0]);
    foreach my $item(@str){
        my $Q=(ord($item)-33);
        $totalitems++;    
        if($Q>$Qfloor)
        {
            $passingitems++;
        }
        
    }
    if($passingitems*100/$percentallowed>=$totalitems&&$totalitems>=$minLength&&$totalitems<=$maxLength)
    {
        return 1;
    }
    else
    {
        return 0;
    }

}
sub getfilename{
	my $str=$_[0];
	my @t=split(/\./,$str);
	return $t[0];	
}
sub getfileext{
	my $str=$_[0];
	my @t=split(/\./,$str);
	return $t[1];	
}
