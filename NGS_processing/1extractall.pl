#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Copy;
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
copyDirectories(1,3,0,$directories[$choice-1]);
closedir $dh;


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
				my $x=3;#used to be 5
				print "+++ ".$dir."/".$_."\n";
				my $file = substr($_,0,4)."_S".substr($_,2,3)."L001_R1_001.fastq.gz";
				print "{{{{{{{{{{{{{".$file."}}}}}}}}}}}\n";
				my $sourcefile=truncDirectory($x,$dir."/".$_."/".$file);
				my $destdir=truncDirectory($x,moveDirectory($origin-$dest,$dir."/".$_)."/Extracted/")."/";
				print "SOURCE ".$sourcefile."\n";
				print "DEST ".$destdir."\n";
				my $sys_call = "7zG x ".$sourcefile." -o".$destdir;
				#my $sys_call = "7zG x ".$sourcefile;
				print $sys_call."\n";
				system($sys_call);#."\n";
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
	#print "---".$output."\t".$delta."\n";
	return substr($output,$delta,length($output)-$delta);

}

#args: name of subfolder (from choice above)
sub unzipAll{
	my $subfolder=shift;
	system("7za a -$subfolder *.fastq");
}
