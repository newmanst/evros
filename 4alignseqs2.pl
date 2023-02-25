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
#
#		edited for PLAseq applications
#
#######################################################################
#               arguments:
#    perl alignseqs.pl
#    ____    this         
#
#  remove sequences longer than 52.
#######################################################################


############## read current directory #####################
my $some_dir=getcwd;
my @directories=();
my @foldernames=();
opendir(my $dh, $some_dir) || die "Can't open $some_dir: $!";
my $c=0;
####line 11###
while (readdir $dh) {
	if(index($_,".")==-1){
		$directories[$c]=$some_dir."/".$_;
		$foldernames[$c]=$_;
		print "".($c+1).") ".$_."\n";
		$c++;
	}
}
print "\nWhich directory would you like to choose?";
my $choice=<STDIN>;
print "\n";
closedir $dh;
###########################################################
my @filenames=getFiles2($some_dir,$foldernames[$choice-1]);
my $in_location="".$foldernames[$choice-1]."/Compiled/";
my $out_location="".$foldernames[$choice-1]."/Aligned/";

############### create new directory if not there ###############
my $new_directory=$foldernames[$choice-1]."/Aligned/";
unless(-e $new_directory or mkdir($new_directory, 0755)) {
        die "Unable to create $new_directory\n";
}

###############################################################################################
################ template selection ###########################################################
###############################################################################################
############################## ENTER ADDITIONAL TEMPLATES HERE ################################
###############################################################################################
#	       1     2     3     4     5     6     7     8     9
my @options=("GFP","TNF","IL1","IL7","ICA","OLD","OLG","RPE","CRP","IL6","ILB","GDF","HF1","HF2");
my %possible_templates=();
$possible_templates{"GFP"}="NNNNNNNNNNGTCTACATCCCTGAGTCTATATGATGTAGATANNNNNNNNNN";#1
$possible_templates{"TNF"}="NNNNNNNNNNTGAGGCCCTATCCCTGAGTCTATATGCCCGATNNNNNNNNNN";#2
$possible_templates{"IL1"}="NNNNNNNNNNTGAGGCCCTATCCCTGAGTCTATATGCCCGATNNNNNNNNNN";#3
$possible_templates{"IL7"}="NNNNNNNNNNAAGTAGTATATGATGACTATCGTACCCAGTAANNNNNNNNNN";#4
$possible_templates{"ICA"}="NNNNNNNNNNCTGCCGTATGATGACTATCGTACCCTGCGACGNNNNNNNNNN";#5
$possible_templates{"OLD"}="GCTATTNNNNNNNNNNTGCCGCGCTTTCGCAGTATTNNNNNNNNNNGATACT";#6
$possible_templates{"OLG"}="NNNNNNNNNNTGTAGAGATGACTATCGTACCCTGCACGGTACNNNNNNNNNN";#7
$possible_templates{"RPE"}="NNNNNNNNNNTCTGGGATGATGACTATCGTACCCTGTAGCACNNNNNNNNNN";#8
$possible_templates{"CRP"}="NNNNNNNNNNTCTGGGATGATGACTATCGTACCCTGTAGCACNNNNNNNNNN";#9
$possible_templates{"IL6"}="NNNNNNNNNNCGTGGGTATCCCTGAGTCTATATGATGACCTCNNNNNNNNNN";#10
$possible_templates{"ILB"}="NNNNNNNNNNAAGTAGTATATGATGACTATCGTACCCAGTAANNNNNNNNNN";#11
$possible_templates{"GDF"}="NNNNNNNNNNCTGCCGTATGATGACTATCGTACCCTGCGACGNNNNNNNNNN";#12
$possible_templates{"HF1"}="NNNNNNNNNNTGTAGAGATGACTATCGTACCCTGCACGGTACNNNNNNNNNN";#13
$possible_templates{"HF2"}="NNNNNNNNNNGTCTACATCCCTGAGTCTATATGATGTAGATANNNNNNNNNN";#14
my %full_name=();
$full_name{"GFP"}="GFP   ";
$full_name{"TNF"}="TNFa  ";
$full_name{"IL1"}="IL-1ra";
$full_name{"IL7"}="IL-7  ";
$full_name{"ICA"}="ICAM-1";
$full_name{"OLD"}="OLD   ";
$full_name{"OLG"}="Oligo ";
$full_name{"RPE"}="R-PE  ";
$full_name{"CRP"}="CRP   ";
$full_name{"IL6"}="IL-6  ";
$full_name{"ILB"}="IL-1B ";
$full_name{"GDF"}="GDF-15";
$full_name{"HF1"}="HiFRe Target ";
$full_name{"HF2"}="HiFRe Control";
my @optionA=(1,2);
my $optionA="GFP TNF";
my @optionB=(1,3,4,5);
my $optionB="GFP IL-1ra IL-7 ICAM-1";
my @optionC=(1,7,8);
my $optionC="GFP R-PE OligoControl";
my @optionD=(1,3,7,9,10,11,12);
my $optionD="GFP OligoControl IL-1ra IL-6 CRP IL-1B GDF-15";
my @optionE=(13,14);
my $optionE="HiFRe Test";
###############################################################################################
###############################################################################################
###############################################################################################

print "\nWhich templates would you like to align to?\n";
print "   "."TARGET\tSEQUENCE\n";
for(my $i=0;$i<scalar(@options);$i++){
    print "".($i+1).") ".$full_name{$options[$i]}."\t".$possible_templates{$options[$i]}."\n";
}
print "please either enter the corresponding numbers separated by spaces \nor choose one of the pre-defined selections bellow:\n";
print "A) ".$optionA."\n";
print "B) ".$optionB."\n";
print "C) ".$optionC."\n";
print "D) ".$optionD."\n";
print "E) ".$optionE."\n\n";
print "ENTER SELECTION: ";
my $options_choice=<STDIN>;
print "\n";
my @targets=();
my $flag=1;
my @IDs=();
my %templates=();
my @result=();
if($options_choice eq "A\n"){
    @result=@optionA;
}elsif($options_choice eq "B\n"){
    @result=@optionB;
}elsif($options_choice eq "C\n"){
    @result=@optionC;
}elsif($options_choice eq "D\n"){
    @result=@optionD;
}elsif($options_choice eq "E\n"){
    @result=@optionE;
}else{
    @result=split(" ",$options_choice);
}
    my $count=0;
    for(my $i=0;$i<scalar(@result);$i++){
        if($result[$i]-1 <= (scalar @options) && exists($possible_templates{$options[$result[$i]-1]})){
            $IDs[$count]=$options[$result[$i]-1];
            $templates{$IDs[$count]}=$possible_templates{$IDs[$count]};
            $count++;
        }else{
            print $result[$i]."is an invalid choice\nWould you like to continue without including ".$result[$i]."?<Y to resume>";
            my $continue=<STDIN>;
            if($continue ne "Y\n"){
                die "goodbye";
            }
        }
    }

print "CHOSEN TEMPLATES:\n";
my $cou=0;
foreach my $k(sort { $a cmp $b} keys %templates){
    #print $full_name{$k}."\t".$IDs[$cou]."\t".$templates{$k}."\n";
    print $full_name{$IDs[$cou]}."\t".$IDs[$cou]."\t".$templates{$IDs[$cou]}."\n";
    $cou++;
}

if(0==0){
open(my $outfile, '>',$foldernames[$choice-1]."/AlignmentOutput.txt");

print "\t\ttotal\tMBC\tmax \tsingle\tdouble\n";
print $outfile "\t\ttotal count\tMBC count\tmax count\tsingle count\tdouble count\n";
    
for(my $q=0;$q<scalar(@filenames);$q++)
{
    #my %outfiles=();
    my @outfiles=();
    my @splitfile=split(/\./,$filenames[$q]);
    
    my %IDread_count=();
    my %IDmbc_count=();
    my %IDmax_count=();
    my %IDsingle_count=();
    my %IDdouble_count=();
    for(my $i=0;$i< scalar @IDs; $i++ ){
        #open($outfiles{$IDs[$i]},'>',$out_location.$splitfile[0]."_".$IDs[$i].".fasta");
        open($outfiles[$i],'>',$out_location.$splitfile[0]."_".$IDs[$i].".fasta");
        $IDread_count{$IDs[$i]}=0;
        $IDmbc_count{$IDs[$i]}=0;
        $IDmax_count{$IDs[$i]}=0;
        $IDsingle_count{$IDs[$i]}=0;
        $IDdouble_count{$IDs[$i]}=0;
    }
    
    my $info="";
    my $seq="";
    my $line_pos=0;
    my $readcount=0;
    open(INPUT, $in_location.$filenames[$q]) || die("\nERROR: Invalid file");
    #print $in_location.$filenames[$q]."\n";
    while(my $in_line=<INPUT>){
        #if($readcount>10000){last;}
        $in_line=~s/\n//g;
        if($line_pos==0){
            $info=$in_line;
            $line_pos++;
        }elsif($line_pos==1){
            $seq=$in_line;
            $readcount++;
            $line_pos=0;
            ###############################################
            #### code the alignment to template here
            #### print to corresponding file!
            ###############################################
            my $MBC1=substr($seq,0,10);
            my $MBC2=substr($seq,-10,10);
            my $max=15;
            my $max_ID_num=-1;
            for(my $i=0;$i< scalar @IDs; $i++ ){
                my $x = similarity($seq,$templates{$IDs[$i]});
                if($x>$max){
                    $max=$x;
                    $max_ID_num=$i;
                }
            }
            if($max_ID_num!=-1){#sequence must align to at least one file
                my @infos=split("-",$info);
                if($infos[1]>=$IDmax_count{$IDs[$max_ID_num]}){
                    $IDmax_count{$IDs[$max_ID_num]}=$infos[1];
                }
                $IDread_count{$IDs[$max_ID_num]}+=$infos[1];
                $IDmbc_count{$IDs[$max_ID_num]}++;
                if($infos[1]==1){
                    $IDsingle_count{$IDs[$max_ID_num]}++;
                }
                if($infos[1]==2){
                    $IDdouble_count{$IDs[$max_ID_num]}++;
                }
                my $fh=$outfiles[$max_ID_num];
                my @temp=split("-",$info);
                my $new_info=">".substr($temp[0],1)."-".$temp[1]."-".$temp[2];
                #my $new_info=">".substr($temp[0],1)."-".$temp[1]."-".$temp[2]."-".substr($seq,10,32);
                print $fh "".$new_info."\n".$MBC1.$MBC2."\n";
		    #print $outfiles{$max_ID} "".$info."\n".$MBC1.$MBC2."\n";
            }
            
        }
    }
    for(my $i=0;$i< scalar @IDs; $i++ ){ # total count, MBC count, max count, single count, double count
        print $filenames[$q]."\t".$IDs[$i]."\t".$IDread_count{$IDs[$i]}."\t".$IDmbc_count{$IDs[$i]}."\t".$IDmax_count{$IDs[$i]}."\t".$IDsingle_count{$IDs[$i]}."\t".$IDdouble_count{$IDs[$i]}."\n";
        print $outfile $filenames[$q]."\t".$IDs[$i]."\t".$IDread_count{$IDs[$i]}."\t".$IDmbc_count{$IDs[$i]}."\t".$IDmax_count{$IDs[$i]}."\t".$IDsingle_count{$IDs[$i]}."\t".$IDdouble_count{$IDs[$i]}."\n";
    }
    
    close(INPUT);
    for(my $i=0;$i< scalar @IDs; $i++ ){
        #close($outfiles{$IDs[$i]});
        close($outfiles[$i]);
    }
}

close($outfile);
}

##################################################################
##################################################################
########################## subroutines ###########################
##################################################################
##################################################################

sub similarity{
    my $one=shift;
    my $two=shift;
    #print "\t\t".$one."\n";
    #print "\t\t".$two."\n";
    my @str1=split(//,$one);
    my @str2=split(//,$two);
    my $c=0;
    for(my $i=10;$i<scalar(@str1)-10;$i++){
        #print "\t\t\t".$str1[$i]."|".$str2[$i]."\n";
        $c += ($str1[$i] eq $str2[$i] ? 1 : 0);
    }
    return $c;
}

#args: counter,start subfolder
sub getFiles2{

	my $d = shift;
	my $f = shift;
	my @fs=();
	my $count=0;
	opendir(my $dh, $d."/".$f."/Compiled") || die "shit";
	while (readdir $dh) {
		if($_ ne "." && $_ ne ".." && $_ ne ".DS_Store"){
			print"\t";
        		print "$_\n";
        		$fs[$count]=$_;
        		$count++;
		}
    	}
	closedir $dh;
	return @fs;
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

