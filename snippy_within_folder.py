#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:30:16 2022
Updated on Mon Oct 3 17:19:20 2022
Updated on Tue Aug 29 16:53:18 2023
@author: nickduggett

Version:1.02


This script will run your samples through snippy, snippy-core and/or snp-dists depending on the option that you provide. It has the option of running on rawreads, assemblies or a mix of these.

Use this line for information:
    python ~/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/Scripts/snippy_within_folder.py -h 

Example command from location of your rawreads and assemblies using E. coli K12 as a reference: 
    python ~/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/Scripts/snippy_within_folder.py -t 3 -r ecoli -m all -f mix
"""

import os, argparse
from argparse import ArgumentParser, SUPPRESS
from pathlib import Path
from shutil import which
from datetime import datetime
import subprocess

###Setting up some basic information as variables that can be called later
today= datetime.now()
date_of_run=today.strftime('%Y%m%d')
output_folder='snippy_'+date_of_run
path='./'
p = Path('~').expanduser()
posixpath_convert= str(p)
fsx_044_path=posixpath_convert+'/mnt/fsx-044/**/*_R1_*fastq.gz'
user=posixpath_convert.split('/')[-1]
path = './'
files = os.listdir(path)
R1_pattern="_R1"
R2_pattern=R1_pattern.replace("R1","R2")
assembly_pattern=[".fa",".fna",".fasta"]

###Adding in arguments that can be called from the commandline. Turning off help
###message so that we can have required and optional arguments in the help that is printed
parser = argparse.ArgumentParser(add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
optional.add_argument('-h','--help',action='help', default=SUPPRESS,help='show this help message and exit')
required.add_argument('--reference','--ref', '-r', help='Path to the reference you want to use. Tip - use "list" to show the available references.', type=str, nargs='+')
optional.add_argument('--threads', '--cpus', '-t', help='How threads you want to use. [Default:1]',type=str, nargs='?', const=1, default="1")
optional.add_argument('--mode', '-m', help='Which mode do you want to run the script in: snippy (runs snippy only), snippy-core (runs snippy and then snippy-core), all (snippy, snippy-core, snp-dists). [Default:snippy]',type=str, nargs='?', const=1, default="snippy")
optional.add_argument('--filetype', '-f', help='What file type do you want to run: fastq (fastqs only), assemblies (assemblies only), mix (mix of fastqs and assemblies). [Default:fastq]',type=str, nargs='?', const=1, default="fastq")
args=parser.parse_args()

###Setting passed arguments as variables
mode=args.mode
threads=args.threads
filetype=args.filetype

refs='\nHere is a list of available references:\n'\
      'Escherichia coli - ecoli, escherichia, Ecoli\n'\
      'Klebsiella pneumoniae - kpneumoniae, klebsiella, Kpneumoniae\n'\
      'Campylobacter jejuni - cjejuni, campylobacter, Cjejuni\n'\
      'Brachyspira hyodysenteriae - bhyo, brachyspirahyodysenteriae, Bhyodysenteriae\n'\
      'Brachyspira intermedia - bint, brachyspiraintermedia, Bintermedia\n'\
      'Brachyspira pilosicoli - bpilo, brachyspirapilosicoli, Bpilosicoli\n'\
      'Brachyspira hampsonii - bhamp, brachyspirahampsonii, Bhampsonii\n'\
      'Brachyspira suanatina - bsuan, brachyspirasuanatina, Bsuanatina\n'
      
###Checking reference has been called
reference=args.reference
reference="".join(map(str,reference))
if args.reference is None:
    print('You need to provide a reference file to the script. See file_search.py -h for details.')
    exit()
elif reference in ['ecoli','escherichia','Ecoli']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/E.coli_MG1655_U00096.3.gbk'
elif args.reference in ['kpneumoniae','klebsiella','Kpneumoniae']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/Klebsiella.pneumoniae.subsp.pneumoniae.HS11286.gbk'
elif args.reference in ['cjejuni','campylobacter','Cjejuni']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/C.jejuni_NCTC11168_ATCC700819_NC002163.gbk'
elif args.reference in ['bint','brachyspiraintermedia','Bintermedia']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/Brachyspira_intermedia_PWS_A_complete_sequence.gbk'
elif args.reference in ['bhyo','brachyspirahyodysenteriae','Bhyodysenteriae']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/Brachyspira_hyodysenteriae_WA1.gbk'    
elif args.reference in ['bpilo','brachyspirapilosicoli','Bpilosicoli']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/Brachyspira_pilosicoli.gbk '    
elif args.reference in ['bhamp','brachyspirahampsonii','Bhampsonii']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/Brachyspira_hampsonii.gbk'    
elif args.reference in ['bsuan','brachyspirasuanatina','Bsuanatina']:
    reference=posixpath_convert+'/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/NCBI_Genomes/Brachyspira_suanatina.gbk'    
elif args.reference == ['list']:
    print(refs)
    exit()
else:
    reference=args.reference
reference="".join(map(str,reference))
print("You have chosen the reference as: "+reference+"\n")
reference_name=reference.split("/")[-1]

    
###Checking mode that has been called
if mode == 'snippy' or mode== '':
    print('You have chosen to run your isolates through fastp and snippy only\n')
elif mode == 'snippy-core':
    print('You have chosen to run your isolates through fastp, snippy and snippy-core\n')
elif mode == 'all':
    print('You have chosen to run your isolates through fastp, snippy, snippy-core and snp-dists\n')
else:
    print(mode+' is not an accepted input for mode\n')
    print('You can only choose a mode from: snippy, snippy-core or all')
    exit()
     
### Checking that the required tools are installed
def is_tool(name):
    """Check whether 'name' is in PATH and marked as executable."""
    return which(name) is not None

tools = ['snippy','fastp','snp-dists']
tool_error='*** Make sure snippy, fastp, snp-dists and multiqc are installed and are in your PATH ***'
tool_install='\nUse this command to install snippy, fastp and snp-dists : mamba install -c bioconda snippy snpEff=5.0 fastp snp-dists multiqc'

if is_tool(tools[0])==False: 
       print('snippy is not installed\n')
       print(tool_error)
       print(tool_install)
       exit()
if is_tool(tools[1])==False:
        print('fastp is not installed\n')
        print(tool_error)
        print(tool_install)
        exit()
if is_tool(tools[2])==False:
        print('snp-dists is not installed\n')
        print(tool_error)
        print(tool_install)
        exit()
else:
        print("All tools are installed, continuing\n")

#Checking that versions of the tools are compatible with the running of the script
def get_version(command):
    version_command = subprocess.check_output(command, stderr=subprocess.STDOUT)
    version = version_command.decode('utf-8').strip().split(" ")[1]
    return int(version.replace(".", ""))

fastp_version = get_version(['fastp', '--version'])
snippy_version = get_version(['snippy', '--version'])

if fastp_version < 234:
    print("Your version of fastp needs to be at least 0.23.4 to run this script\n\nPlease install the latest version before running")
    exit()

if snippy_version < 460:
    print("Your version of snippy needs to be at least 4.6.0 to run this script\n\nPlease install the latest version before running")
    exit()
 
#Determining what the input files are and dealing with them as required       
if filetype == 'fastq' or filetype == '':
    with open(path+"sample_list.txt","w") as sample_files:
        for fastq_R1 in files:
            if fastq_R1.__contains__(R1_pattern)&fastq_R1.endswith("fastq.gz"):
                fastq_R2=fastq_R1.replace(R1_pattern,R2_pattern)
                sample_files.write(fastq_R1+"\t"+fastq_R2+'\n')
elif filetype == 'assemblies':
    with open(path+"sample_list.txt", "w") as sample_files:
        for assembly in files:
            if assembly.endswith(tuple(assembly_pattern)):
                sample_files.write(assembly + '\n')
elif filetype == 'mix':
    with open(path+"sample_list.txt","w") as sample_files:
        for fastq_R1 in files:
            if fastq_R1.__contains__(R1_pattern)&fastq_R1.endswith("fastq.gz"):
                fastq_R2=fastq_R1.replace(R1_pattern,R2_pattern)
                sample_files.write(fastq_R1+"\t"+fastq_R2+'\n')
    with open("sample_list.txt", "a") as sample_files:
        for assembly in files:
            if assembly.endswith(tuple(assembly_pattern)):
                sample_files.write(assembly + '\n')
else:
    print(filetype+' is not an accepted input for filetype\n')
    print('You can only choose a filetype from: fastq, assemblies or mix')
    exit()
           
list_of_samples = open("sample_list.txt","r")
lines=list_of_samples.readlines()

# Check if the file is empty
if len(lines) == 0:
    print("I can't find any files to process. Exiting...")
    list_of_samples.close()
    exit()

#Process each line of the sample file which equates to one line per sample
for line in lines:
    R1=line.split("\t")[0]
    if R1.endswith(".fastq.gz"):
        R2=line.split("\t")[1]
        R2=R2.replace("\n","")
        fastp_filename_R1=R1.replace(".fastq","_filtered.fastq")
        fastp_filename_R2=R2.replace(".fastq","_filtered.fastq")
    else:
        R1=R1.replace("\n","")
    isolate_name=R1.split("_")[0]
    snippy_output_name=R1.split("_")[0]+"_snippy"
    snippy_output_name_assembly=R1.split(".")[0]+"_snippy"
    checkFile=snippy_output_name+'/snps.vcf'
    checkFile_assembly=snippy_output_name_assembly+'/snps.vcf'
    #If the sample is a raw read then process it like this
    if R1.endswith(".fastq.gz"):
        run_fastp="fastp -i "+R1+" -I "+R2+" -o "+fastp_filename_R1+" -O "+fastp_filename_R2+" -l 80 -r --cut_right_window_size 10 -w "+threads+" -j /dev/null -h /dev/null"
        run_snippy="snippy --R1 "+fastp_filename_R1+" --R2 "+fastp_filename_R2+" --ref "+reference+" --cpus "+threads+" --outdir "+snippy_output_name+" --force  --cleanup"
        cleanup="rm "+fastp_filename_R1+" "+fastp_filename_R2
        print("\nRunning fastp on "+R1.split("_")[0]+" with "+threads+" threads")
        os.system(run_fastp)
        print("snippy --R1 "+fastp_filename_R1+" --R2 "+fastp_filename_R2+" --ref "+reference+" --cpus "+threads+" --outdir "+snippy_output_name+" --force  --cleanup")
        print("\nRunning snippy on "+R1.split("_")[0]+" with "+threads+" threads")
        os.system(run_snippy)
        if os.path.isfile(checkFile):
            print("\nSnippy worked for "+R1.split("_")[0])
        else:
            print("\nSnippy failed to produce the snps.vcf file, check the snps.log file for information\n\nExiting")
            exit()
        print("\nCleaning up files generated by fastp")
        os.system(cleanup)
    #Else if the sample is an assembly, process it like this    
    else:
        run_snippy="snippy --ctgs "+R1+" --ref "+reference+" --cpus "+threads+" --outdir "+snippy_output_name_assembly+" --force  --cleanup"
        print("\nRunning snippy on "+R1.split(".")[0]+" with "+threads+" threads")
        os.system(run_snippy)
        if os.path.isfile(checkFile_assembly):
            print("\nSnippy worked for "+R1.split(".")[0])
        else:
            print("\nSnippy failed to produce the snps.vcf file, check the snps.log file for information\n\nExiting")
            exit()
            
os.system("rm sample_list.txt")
print("\nSummary output of fastp that was used by snippy for mapping can be found at: your_isolate_name.fastp.html")
print("\nTo run snippy-core use the following command: \nsnippy-core --ref "+reference+" *_snippy")

end_message='\nYay, you ('+user+') ran the code and I (the computer), did all the work!'
run_snippy_core='snippy-core --ref '+reference+' --prefix '+date_of_run+'_core '+path+'/*_snippy' #&& mv '+date_of_run+'_core*aln '+path+' && mv '+date_of_run+'_core.* '+path+' '+date_of_run+'*_core*'
run_snp_dists_full_align='snp-dists -j '+threads+' '+path+'/'+date_of_run+'_core.full.aln > '+path+'/snp_dists_'+date_of_run+'_core.full.aln.tsv'
run_snp_dists_core_align='snp-dists -j '+threads+' '+path+'/'+date_of_run+'_core.aln > '+path+'/snp_dists_'+date_of_run+'_core.aln.tsv'
if mode=='snippy':
    print(end_message)
    exit()
elif mode=='snippy-core':
    print('\nRunning snippy-core\n')
    os.system(run_snippy_core)
    print('\nRan snippy-core')
    print(end_message)
    exit()
elif mode=='all':
    print('\nRunning snippy-core\n')
    os.system(run_snippy_core)
    print('\nRunning snp-dists on core.full.aln \n')
    os.system(run_snp_dists_full_align)
    print('\nRunning snp-dists on core.aln \n')
    os.system(run_snp_dists_core_align)
    print(end_message)
    print("\nTo open the whole genome alignment snp-dists file in libreoffice use: libreoffice --calc snp_dists_"+date_of_run+"_core.full.aln.tsv\n")
    print("Or to open the core genome alignment snp-dists file in libreoffice use: libreoffice --calc snp_dists_"+date_of_run+"_core.aln.tsv\n")
    exit()
