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
use Time::HiRes qw( gettimeofday);
my $delimiter="_";
my $help="Help! Something has gone wrong";

#######################################################################
#               arguments:
#    perl compileseqs.pl        
#    ____    this         
#
#    collapses exact copies into a single read and perserves the read count
#######################################################################





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
my @filenames=getFiles2($some_dir,$foldernames[$choice-1]);
my $in_location="".$foldernames[$choice-1]."/Filtered/";
my $out_location="".$foldernames[$choice-1]."/Compiled/";

############### create new directory ###############
my $new_directory=$foldernames[$choice-1]."/Compiled/";
unless(-e $new_directory or mkdir($new_directory, 0755)) {
        die "Unable to create $new_directory\n";
}

for(my $q=0;$q<scalar(@filenames);$q++)
{
    my @temp=split("_",$filenames[$q]);
    my $newfilename = $temp[0].".fasta";
	my $c = "perl SupportFiles/fastaptamer_count.pl -i \"".$in_location.$filenames[$q]."\" -o \"".$out_location.$newfilename."\"";
	print "".$c."\n";
	system($c);
}



##################################################################
##################################################################
########################## subroutines ###########################
##################################################################
##################################################################


#args: counter,start subfolder
sub getFiles2{

	my $dir = shift;
	my $folder = shift;
	my @files=();
	my $count=0;
	opendir(my $dh, $dir."/".$folder."/Filtered") || die "oops";
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
