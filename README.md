# EVROS Data Processing 
Source code for the work to "Extending the dynamic range of biomarker quantification through molecular equalization" presented here: https://www.biorxiv.org/content/10.1101/2023.02.07.527534v1

The files here do:
1) **qPCR_processing**: The qPCR raw processing to Cqs and present demonstration of converting to DNA concentrations given a calibration file 
2) **NGS_processing**: The sequence parsing of the fastq files 

The steps for NGS processing are as follows: 

## STEP 1: Folder set up. 
Establish folder that will contain all sequencing data. I usually save them in experiment folder with a name “NGS”.  Copy perl files into folder. Make folder for each sequencing experiment. Download fastqs from Basespace. Put raw sequencing files into those folders. (Note you can update the perl script to not have to be copied this way) 

## STEP 2: Extract and unzip all fastq.gz files . 
This step is only for those on MacOS. This can be done automatically with 1extractall.pl, if on Windows. Extracted files must be put into folder titled “Extracted”.

2a. Go to the folder where just before all folder names “A1_L001...”, “A2_L001” ...  (probaby labeled as “FASTQ_Generation_<DATE> .... “. Use commands such as (cd means change directory)
    
    cd <top-directory>/NGS/<directory>  

2b. move all files in the subdirectory with R1 reads to the current one 
 
    find . -name "*R1_001.fastq.gz" -exec mv -i -- {} . \;

2c. unzip all these files: (code unzips all folders with gz as a end)

    gunzip *.gz	 

Note: If you have output in terminal such as "<...>.fastq.qz has 1 other links -- skipping". Force it to unzip by running: 
    
    gunzip *.gz -f 
  
2d. Remove all files in other folder replace <FASTQ_Generation_folder> with actual file name: (*ATTENTION*: this removes all folders remaining in the called folder!!!! Make sure you have a backup)
    
    rm -rf <FASTQ_Generation_folder>
  
2e. Because perl files require a specific folder naming hierarchy of the fastq files, we put all the extracted files into a file named “Extracted” 
  
    mkdir Extracted 
    mv *.fastq Extracted/   
 
## STEP 3: Filter reads 

3a. change directory one level out to the folder with your perl scripts 

    cd ../
  
3b. Run “perl 2filterseqs.pl min_length max_length min_percent min_Q” 
    min_length/max_length — apply limits to sequence length (default 52,52)
   -  Only considers sequences with length between min_length and max_length and quality where at least min_percent% of bases have Q>=min_Q.
   - Default args: 52 52 100 20
   - Choose your experiment. Output will be stored in “Filtered” folder
   - If default arguments:
  
  perl 2filterseqstrim9.pl
  
## STEP 4: compile sequences into unique reads and count.
Run “perl 3compileseqs.pl”. Choose your experiment. Output will be stored in “Compiled” folder

*note I added here “SupportFiles” folder which holds a fastaptamer_count.pl folder required to run 4. Ideally this would already be set up earlier in step 1 when importing code.  
  
## STEP 5: Align sequences 
Align sequences to known reporter sequences to determine which reporter generated that sequence. *This step takes the longest esp when run locally on your computer.* 
   
 Run “perl 4alignseqs.pl”. Choose your experiment. Choose the template desired Output will be stored in “Aligned” folder.
 The .txt output titled “AlignmentOutput.txt” will be output in experiment folder.
 Note: You will need to update the sequences used if you have different sequences for your probes here. 
  
Congrats! Your NGS data is now pre-processed to use for counting! 
  
  

