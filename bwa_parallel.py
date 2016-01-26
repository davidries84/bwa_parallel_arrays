# 0. einlesen, Variablen definieren, Directories erstellen

import os
import sys
import glob
import commands
import socket
import smtplib
from email.mime.text import MIMEText
import re
from optparse import OptionParser
import shutil
from subprocess import call, Popen, PIPE
import time
import datetime 
from itertools import izip_longest

def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def calcSplitLen(file, blocks, chunk = 4):
        '''calculates the number of lines per block, for splitting file in number of blocks files of approx. eqal size. chunk is the number of lines that must not be split, for fastq files this is 4 = 4'''
	# calc number of lines
	print "Calculating file length"
        cmd = "wc -l " + file # fastest way
        p = Popen(cmd, stdout=PIPE, shell=True)
        p.wait()
        out, err = p.communicate()
        filelength = int(out.split(" ")[0])
        div = filelength / split_number # div =number of lines per split file
        print "filelength: " + str(filelength)
        if div < chunk:
            splitLen = chunk
        elif div % chunk == 0:
            splitLen = div
        else:
            #mod = div % chunk
            #add = chunk - mod
            #splitLen = div + add
            splitLen = div + chunk
          
        return splitLen

def split_file(file,splitLen,split_dir):
        '''splits the file in files of length splitLen'''	
	outfile = os.path.basename(file).replace('.fastq','_split_{0}.fastq')
        outfiles = []
       
	with open(file) as f:
    		for i, g in enumerate(grouper(splitLen, f, fillvalue=''), 1):
    			outfiles.append(outfile.format(i))
        		with open(outfile.format(i), 'w') as fout:
        			# B2444_ATCACG_L008_R1_001.paired.fastq
        			# B2444_ATCACG_L008_R1_001.paired_split_40.fastq
            			fout.writelines(g)
	return outfiles

def time_difference(sec):
    if sec < 60:
	    return "0 day(s), 0 hour(s), 0 minute(s), " + str(sec) + " second(s) [" + str(sec) + " seconds]"
    else:
        days = int(sec / 86400)
        hours = int((sec / 3600) - (days * 24))
        minutes = int((sec / 60) - (days * 1440) - (hours * 60))
        seconds = int(sec % 60) 
        return str(days) + " day(s), " + str(hours) + " hour(s), " + str(minutes) + " minute(s), " + str(seconds) + " second(s) [" + str(sec) + " seconds]"


def write_array_file(path_to_bwa, path_to_samtools, path_to_reference, split_file, path_to_splitFiles, threads, splits):
	'''create and write the array file containing the commands to map fastq files to a reference creating bams'''


	filename = os.path.basename(split_file)
	fw_splitted_files = os.path.join(path_to_splitFiles,filename.replace('.fastq','_split_$1.fastq'))
	rv_splitted_files = fw_splitted_files.replace('R1','R2')
	bamFile = os.path.join(os.path.dirname(fw_splitted_files) ,os.path.basename(fw_splitted_files).replace('.fastq', '.bam'))
	mapping_cmd = path_to_bwa + " mem  -t " + str(threads) + " " + path_to_reference + " " + os.path.abspath(fw_splitted_files) + " " + os.path.abspath(rv_splitted_files) + " | " + "/vol/biotools/bin/samtools view -bS - |  /vol/biotools/bin/samtools sort -m 200000000 - " +  bamFile.replace('.bam','')

	#writing shell script
	filehandle = open('mapping.sh','w')
	os.chmod('mapping.sh', 0700)
	filehandle.writelines('#!/bin/bash \n')
	filehandle.writelines(mapping_cmd)
	filehandle.close()

	#writing array, that runs shell script N times

	array_cmd = os.path.abspath('mapping.sh')  + ' $SGE_TASK_ID'
	filehandle = open('mapping_array.csh', 'w')
	os.chmod('mapping_array.csh', 0700)
	filehandle.writelines(array_cmd)
	filehandle.close()
        

	# run array
	# cmd = qsub -t 1-n -N mapping -vf 8G -l arch=lx24-amd64 -cwd -pe multislot "Threads"-e " + error_log + " -o output mapping_array.csh
	#calculate memory requirements
	virtualMem = int(8)/int(threads)
	cmd = "qsub -t 1-" + str(int(splits)) + " -N mapping -l idle=1 -l arch=lx24-amd64,vf=" + str(virtualMem) + "G -cwd -sync y -pe multislot " + str(threads) + " -e  error_log  mapping_array.csh"
        p = Popen(cmd, shell=True, stdin=PIPE,stdout=PIPE)
        array_job_ID = p.communicate()[0].split()[2] # save the job ID 
	print "submitted " 
	print cmd
	print "with Job ID: " + 	array_job_ID
#	p.wait() # wait for the array to finish
	path_to_mapped_files = os.path.dirname(bamFile)
	
	return array_job_ID.split('.')[0], path_to_mapped_files, p

def write_array_file_singlereads(path_to_bwa, path_to_samtools, path_to_reference, name_split_file, path_to_splitFiles, threads, number_splits):
	'''create and write the array file containing the commands to map fastq files to a reference creating bams. split_file only defines the filename.'''


	filename = os.path.basename(name_split_file)
	fw_splitted_files = os.path.join(path_to_splitFiles,filename.replace('.fastq','_split_$1.fastq'))
	#rv_splitted_files = fw_splitted_files.replace('R1','R2')
	bamFile1 = os.path.join(os.path.dirname(fw_splitted_files) ,os.path.basename(fw_splitted_files).replace('.fastq', '.bam'))
	mapping_cmd1 = path_to_bwa + " mem  -t " + str(threads) + " " + path_to_reference + " " + os.path.abspath(fw_splitted_files) + " | " + "/vol/biotools/bin/samtools view -bS - |  /vol/biotools/bin/samtools sort -m 200000000 - " +  bamFile1.replace('.bam','')
	#bamFile2 = os.path.join(os.path.dirname(rv_splitted_files) ,os.path.basename(rv_splitted_files).replace('.fastq', '.bam'))
	#mapping_cmd2 = path_to_bwa + " mem  -t " + str(threads) + " " + path_to_reference + " " + os.path.abspath(rv_splitted_files) + " | " + "/vol/biotools/bin/samtools view -bS - |  /vol/biotools/bin/samtools sort -m 200000000 - " +  bamFile2.replace('.bam','')


 
	#writing shell script
	filehandle = open('mapping.sh','w')
	os.chmod('mapping.sh', 0700)
	filehandle.writelines('#!/bin/bash \n')
	filehandle.writelines(mapping_cmd1)
	#filehandle.writelines(mapping_cmd2)
	filehandle.close()

	#writing array, that runs shell script N times

	array_cmd = os.path.abspath('mapping.sh')  + ' $SGE_TASK_ID'
	filehandle = open('mapping_array.csh', 'w')
	os.chmod('mapping_array.csh', 0700)
	filehandle.writelines(array_cmd)
	filehandle.close()
        

	# run array
	# cmd = qsub -t 1-n -N mapping -vf 8G -l arch=lx24-amd64 -cwd -pe multislot "Threads"-e " + error_log + " -o output mapping_array.csh
	#calculate memory requirements
	virtualMem = int(8)/int(threads)
	cmd = "qsub -t 1-" + str(int(number_splits)) + " -N mapping -l idle=1 -l arch=lx24-amd64,vf=" + str(virtualMem) + "G -cwd -sync y -pe multislot " + str(threads) + " -e  error_log  mapping_array.csh"
        p = Popen(cmd, shell=True, stdin=PIPE,stdout=PIPE)
        array_job_ID = p.communicate()[0].split()[2] # save the job ID 
	print "submitted " 
	print cmd
	print "with Job ID: " + 	array_job_ID
#	p.wait() # wait for the array to finish
	path_to_mapped_files = os.path.dirname(bamFile1)
	
	return array_job_ID.split('.')[0], path_to_mapped_files, p
	
	
def mergeBams(path_to_samtools, split_file, mappedFilesList ,array_job_ID, mapping_process):
    """merges all bam files MappedFilesList into one"""
    if "unpaired" in split_file:
        final_file = split_file.replace('_001.unpaired.fastq','_unpaired.bam')
        
    else:
        final_file = split_file.replace('_R1_001.paired.fastq','_paired.bam')
    if len(mappedFilesList) == 0:
            print "no mapped files ->  exiting"
            sys.exit()

    elif len(mappedFilesList) == 1:
            print 'waiting for mapping to finish'
            mapping_process.wait()
            cmd = "cp " + mappedFilesList[0] + " " + final_file
            print cmd
            print "merge not needed. waiting for mapping to finish and copying mapped file"
            p = Popen(cmd, shell=True)
            p.wait()
    else:
            files = ''
            files += ' '.join(mappedFilesList) + ' '
            cmd = "qsub -cwd -N bwa_merge -hold_jid  " +  array_job_ID + " -sync y -e " + error_log + " -o " + out_log + " -b y  -l vf=16G " + path_to_samtools + " merge " + final_file+ " " + files
            if not os.path.isfile(final_file):
                    print "merging the sorted bams to " + final_file
                    print cmd
                    p = Popen(cmd, shell=True)
                    p.wait()
            else: 
                    print "samtools merge output already exists"
    
    return final_file

def splitfiles(forward_file, reverse_file, split_size, split_number, split_dir, checks = True, chunk = 4):
    """splits fastq files in mappable chunks. for unpaired reads, leave reverse_file empty."""
    os.chdir(split_dir) # cd to the directory where the split files will be created
    splitfiles = []
    unpaired = False
    if reverse_file == []:
            unpaired=True
            
    if unpaired == False:

            chunk = 4 # chunk length
            splitfiles = []
            files_already_split = False
### calculate split length/ number of lines for each split file
            if split_size > 0 and split_number > 1: 
                print "You have to choose between -s and -n option"
                sys.exit()
            elif split_size > 0:
                div = split_size # lines per file
            elif split_number == 1: ## if no split, create a symbolic link, to save the wc and the copying, but keep compatibility
                forw_sym = os.path.join(os.getcwd(),os.path.basename(forward_file.replace('.fastq','_split_1.fastq')))
                rev_sym = os.path.join(os.getcwd(),os.path.basename(reverse_file.replace('.fastq','_split_1.fastq')))

                print 'Split number is 1. Setting symbolic links instead of splitting.'
                print forw_sym
                print rev_sym


                if not os.path.islink(forw_sym):
                    os.symlink(forward_file, forw_sym)
                else:
                    print "link to forward file already exists, using old sym_link"

                if not os.path.islink(rev_sym):
                    os.symlink(reverse_file, rev_sym)
                else:
                    print "link to reverse file already exists, using old sym_link"

                splitfiles=[forw_sym,rev_sym]
                return splitfiles
            elif split_number == len(glob.glob(forward_file.replace('.fastq',"*split*.fastq"))) + len(glob.glob(reverse_file.replace('.fastq',"*split*.fastq"))) / 2:
                files_already_split = True
            elif split_number > 0:
                files_already_split = False
                div = calcSplitLen(forward_file, split_number)
# div =number of lines per split file
              
                print "lines per split file: " + str(div)

            else:
                print "You have to specify -s or -n (positive int value)"
                sys.exit()

    if unpaired == True:

            files_already_split = False

            if split_size > 0 and split_number > 1: 
                print "You have to choose between -s and -n option"
                sys.exit()
            elif split_size > 0:
                div = split_size # lines per file
            elif split_number == 1: ## if no split, create a symbolic link, to save the wc and the copying, but keep compatibility
                forw_sym = os.path.join(os.getcwd(),os.path.basename(forward_file.replace('.fastq','_split_1.fastq')))
                print 'Split number is 1. Setting symbolic link instead of splitting.'
                print forw_sym
                if not os.path.islink(forw_sym):
                    os.symlink(forward_file, forw_sym)
                else:
                    print "link already exists, using old sym_link"
                
                splitfiles.append(forw_sym)
                return splitfiles

            elif split_number == len(glob.glob(forward_file.replace('.fastq',"*split*.fastq"))):
                files_already_split = True
            elif split_number > 0:
                files_already_split = False
                div = calcSplitLen(forward_file, split_number)
               
                print "lines per split file: " + str(div)
            else:
                print "You have to specify -s or -n (positive int value)"
                sys.exit()

	# splitting
    # check optional (handle with care) not implemented: 

    if files_already_split == False and unpaired == False:
	splitforwfiles = split_file(forward_file,div,split_dir)
        print "split of paired file " + forward_file + " finished"
      	splitrevfiles = split_file(reverse_file,div,split_dir)
        print "split of paired file " + reverse_file + " finished"
        for i, (a,b) in enumerate( zip(splitforwfiles, splitrevfiles)):
            splitfiles.append(a)
            splitfiles.append(b)

    elif files_already_split == False  and unpaired == True:
	splitfiles = split_file(forward_file,div,split_dir)
        print "split of unpaired file " + forward_file + " finished"
       
        
    elif files_already_split == True and unpaired == False :
	splitfiles = zip(glob.glob(forward_file.replace('.fastq',"*split*.fastq")), glob.glob(reverse_file.replace('.fastq',"*split*.fastq")))
	splitfiles = [j for i in splitfiles for j in i]
	print "Found already split files:"
	print splitfiles
    elif files_already_split == True and unpaired == True :
	splitfiles = glob.glob(forward_file.replace('.fastq',"*split*.fastq"))

	print "Found already split files:"
	print splitfiles

	
    else:
	print "split failed"
	    
    os.chdir(os.pardir)
    return splitfiles # immer abwechselnd ein forward und ein reverse split file




""" This script starts bwa mem for all Fastq files in the folder given in base_path. 
The Files have to be named in standard GA2x or Hiseq format i.e. Name_NoIndex_L001_R1_001.fastq. Splitting the files
seems to be the bottleneck."""


init_start = time.time()
# Parsing the command line options
parser = OptionParser(usage="usage: python %prog [options]")
parser.add_option('-b', '--bwapath', dest="bwaPath", type="string", default="/vol/biotools/bin/bwa0_7_5", help='Path to the bwa-Tool. Bwa version must support bwa mem <reference> <in1.fq> <in2.fq>. Default is /vol/biotools/bin/bwa0_7_5.')
parser.add_option('-P', '--basePath', dest="basePath", type="string", default="./", help='Folder containing the Temp-Data')
parser.add_option('-I', '--ioPath', dest="ioPath", type="string", default="./", help='Folder containing the HiSeq raw-Data and the result')
parser.add_option('-F', '--pairedfilter', dest="pairedfilter", type="string", default='*R1*.paired.fastq', help='Filter for paired-end reads raw-data files by filename, if you do not want all the files contained in the folder specified in -P, allows regular expressions')
parser.add_option('-f', '--unpairedfilter', dest="unpairedfilter", type="string", default='*.unpaired.fastq', help='Filter for unpaired reads raw-data files by filename, if you do not want all the files contained in the folder specified in -P, allows regular expressions')
parser.add_option('-R', '--reference', dest="reference", type="string", default="./", help='the reference as fasta file')
parser.add_option('-t', '--mappingthreads', dest="mappingthreads", type="int", default=1, help='number of threads bwa should use')
parser.add_option('-c', action="store_true", dest="checks", default=False, help="enable fastq files correctly paired check ")
parser.add_option('-e', '--errorLog', dest="errorLog", type="string", default="./", help='Folder containing the cluster error log')
parser.add_option('-o', '--outLog', dest="outLog", type="string", default="./", help='Folder containing the cluster output log')
parser.add_option('-s', '--splitsize', dest="splitsize", type="int", help='size of each split, recommended over -n for speed')
parser.add_option('-n', '--splitnumber', dest="splitnumber", type="int", default=1, help='number of splits')
parser.add_option('-T', action="store_true", dest="test", default=False, help="enter testing mode for evaluation (handle with care)")

(options, args) = parser.parse_args()

# the bwa path
bwa_path = os.path.abspath(options.bwaPath)
print "Bwa path: " + bwa_path

# the basepath of the temp data files
base_path = os.path.abspath(options.basePath)
print "Temp path: " + base_path

# the basepath of the raw data files and result
io_path = os.path.abspath(options.ioPath)
print "IO path: " + io_path

# the filename filter pattern
pairedfilter = options.pairedfilter
unpairedfilter = options.unpairedfilter
print "paired-end filename filter: " + pairedfilter
print "unpaired filename filter: " + unpairedfilter

# the reference path
reference = os.path.abspath(options.reference)
print "reference path: " + reference

# the cluster error log
error_log = os.path.abspath(options.errorLog)
print "Cluster error log: " + error_log

# the cluster output log
out_log = os.path.abspath(options.outLog)
print "Cluster output log: " + out_log

# mapping parameters
mapping_threads = str(options.mappingthreads)
print "mapping threads: " + str(options.mappingthreads)

# checks
checks = options.checks
print "Checks: " + str(checks)

# testen
testing = options.test
print "Testing: " + str(testing)

# split size
split_size = options.splitsize
print "split size: " + str(options.splitsize)

# split number
split_number = options.splitnumber
print "number of splits: " + str(options.splitnumber)

forward_files=[]
# list of the forward reads containing fastq files to process with trimmomatic
forward_files=glob.glob(os.path.join(io_path, pairedfilter))

# check for already existing bam files
for filename in reversed(forward_files):
    if os.path.exists(filename.replace('_R1_001.paired.fastq','_paired.bam')):
        print "mapped file for " + filename + " already found"
        forward_files.remove(filename)


if not forward_files:
    print "List is empty. No paired end input files found"

# create list of corresponding reverse files
reverse_files=[]

if forward_files:
        for filename in forward_files:
            if  os.path.exists(filename.replace('R1','R2')):
                print 'Found reverse read file: ' + filename.replace('R1','R2')
                reverse_files.append(filename.replace('R1','R2'))
            else:
                print filename.replace('R1','R2') + ' does not exist'
                print "Reverse files not found. exiting"
                sys.exit()
        forward_files.sort()
        reverse_files.sort()




#create list of unpaired files

unpaired_files=[]
# list of the forward reads containing fastq files to process with trimmomatic
unpaired_files=glob.glob(os.path.join(io_path, unpairedfilter))
print "Unpaired files found: " 
print unpaired_files

# check for already existing bam files
for filename in reversed(unpaired_files):
    if os.path.exists(filename.replace('_001.unpaired.fastq','_unpaired.bam')):
        print "mapped file for " + filename + " already found"
        unpaired_files.remove(filename)


if not unpaired_files:
   print "List is empty. No unpaired files found"

# abort if no files found
if not forward_files and not unpaired_files:
        print "No paired or unpaired files found. exiting"
        sys.exit()

# create lists of all input files
input_files = forward_files + reverse_files

# make index if not made
if not os.path.isfile(reference+'.bwt'):
    # bwa will make the index in the correct dir
    print "Indexing the reference"
    cmd = [bwa_path, 'index', reference]
    call(cmd) # generate index and wait until it's done
    
if not os.path.isfile(reference+'.ann'):
    # bwa will make the index in the correct dir
    print "Indexing the reference for bwtsw"
    cmd = [bwa_path, 'index', '-a', 'bwtsw', reference]
    call(cmd) # generate index and wait until it's done



# Directory fuer cluster output log erstellen
if not os.path.exists(out_log):
    os.mkdir(out_log)

# Directory fuer cluster error log erstellen
if not os.path.exists(error_log):
    os.mkdir(error_log)

init_end = time.time()


###mapping unpaired data
for i in range(0, len(unpaired_files)):

    print "processing: "
    print unpaired_files[i]

    # Directory fuer temp files erstellen
    if not os.path.exists(base_path):
        os.mkdir(base_path, 0755)

    # Directory fuer splitted files erstellen
    if not os.path.exists(base_path+"/splitted_files"):
        os.mkdir(base_path+"/splitted_files", 0755)

    # Directory fuer sam files erstellen
    if not os.path.exists(base_path+"/sam_files"):
        os.mkdir(base_path+"/sam_files", 0755)

    # Directory fuer mapped files erstellen
    if not os.path.exists(base_path+"/mapped_files"):
        os.mkdir(base_path+"/mapped_files", 0755)

    # Directory fuer sorted mapped files erstellen
    if not os.path.exists(base_path+"/sorted_mapped_files"):
        os.mkdir(base_path+"/sorted_mapped_files", 0755)

    splitten_start = time.time()
    # 1. splitten
    print "splitting. this may take some time... "
    splitfilesList = splitfiles(unpaired_files[i], [], split_size, split_number,base_path+"/splitted_files/" , checks)
    print 'split files to:'
    print splitfilesList
    splitten_end = time.time()
    print "Time for split: " + str(time_difference(splitten_end - splitten_start))

    # 2. mappen
    # Listen der splitted files
    print "mapping... "
    
    array_job_ID, path_to_mapped_files, mapping_process = write_array_file_singlereads(bwa_path, "/vol/biotools/bin/samtools", reference, unpaired_files[i], base_path+"/splitted_files", mapping_threads, len(splitfilesList))
    mappen_end = time.time()

    print 'files mapped to:'
    print path_to_mapped_files

    # 3. mergen
    
    mappedFilesList =  []
    for j in range(0,len(splitfilesList)):
	    mappedFilesList.append(os.path.join(path_to_mapped_files , splitfilesList[j].replace('.fastq','.bam')))

    print 'files mapped to:'
    print mappedFilesList

	    
    merged_bam = mergeBams("/vol/biotools/bin/samtools", unpaired_files[i] , mappedFilesList ,array_job_ID, mapping_process)
    print "merged bams to : " + merged_bam
    mergen_end = time.time()


    # 4. cleanup
    cleanup = False
    def purge(dir, pattern):
        for f in os.listdir(dir):
            if re.search(pattern, f):
                os.remove(os.path.join(dir, f))

    os.chdir(io_path)
    if testing:
        purge(io_path, os.path.split(final_file))
        if os.path.exists(base_path):
              print base_path + " exists"
              shutil.rmtree(base_path)
              print "removed"
              purge(error_log, "bwa_mapping*")
              purge(error_log, "bwa_merge*")
              purge(out_log, "bwa_mapping*")
              purge(out_log, "bwa_merge*")
    else:
                if os.path.exists(base_path+"/splitted_files") and cleanup == True:
                    print base_path+"/splitted_files exists"
                    shutil.rmtree(base_path+"/splitted_files")
                    print "removed"
                if os.path.exists(base_path+"/sam_files"):
                    print base_path+"/sam_files exists"
                    shutil.rmtree(base_path+"/sam_files")
                    print "removed"
                if os.path.exists(base_path+"/mapped_files"):
                    print base_path+"/mapped_files exists"
                    shutil.rmtree(base_path+"/mapped_files")
                    print "removed"
                if os.path.exists(base_path+"/sorted_mapped_files"):
                    print base_path+"/sorted_mapped_files exists"
                    shutil.rmtree(base_path+"/sorted_mapped_files")
                    print "removed"


    cleanup_end = time.time()


    dest = open(os.path.join(io_path, os.path.split(unpaired_files[i])[1].replace('_001.unpaired.fastq', "_" + str(split_size) + "s_" + str(split_number) + "n_" + str(mapping_threads) + "t.unpaired.txt")), 'w')

    stats = "Stats to " + unpaired_files[i] + ":\n"
    ss = "Split size: " + str(split_size) + "\n"
    sn = "Number of splits: " + str(split_number) + "\n"
    mt = "Mapping threads: " + str(mapping_threads) + "\n"
    ch = "Checks: " + str(checks) + "\n"
    init = "Init duration: " + time_difference(init_end - init_start) + "\n"
    splitten = "Splitting duration: " + time_difference(splitten_end - splitten_start) + "\n"
    mappen = "Mapping duration: " + time_difference(mappen_end - splitten_end) + "\n"
    mergen = "Merging duration: " + time_difference(mergen_end - mappen_end) + "\n"
    cleanup = "Clean up duration: " + time_difference(cleanup_end - mergen_end) + "\n"
    total = "Total duration: " + time_difference(cleanup_end - init_start) + "\n"

    dest.write(stats)
    dest.write(ss)
    dest.write(sn)
    dest.write(mt)
    dest.write(ch)
    dest.write(init)
    dest.write(splitten)
    dest.write(mappen)
    dest.write(mergen)
    dest.write(cleanup)
    dest.write(total)

    dest.close()




### mapping paired-end files

for i in range(0, len(forward_files)):
    init_start = time.time()
    print "processing: "
    print forward_files[i]
    print reverse_files[i]

    # Directory fuer temp files erstellen
    if not os.path.exists(base_path):
        os.mkdir(base_path, 0755)

    # Directory fuer splitted files erstellen
    if not os.path.exists(base_path+"/splitted_files"):
        os.mkdir(base_path+"/splitted_files", 0755)


    # Directory fuer sam files erstellen
    if not os.path.exists(base_path+"/sam_files"):
        os.mkdir(base_path+"/sam_files", 0755)

    # Directory fuer mapped files erstellen
    if not os.path.exists(base_path+"/mapped_files"):
        os.mkdir(base_path+"/mapped_files", 0755)

    # Directory fuer sorted mapped files erstellen
    if not os.path.exists(base_path+"/sorted_mapped_files"):
        os.mkdir(base_path+"/sorted_mapped_files", 0755)

    splitten_start = time.time()
    # 1. splitten
    print "splitting. this may take some time... "
    splitfilesList = splitfiles(forward_files[i], reverse_files[i], split_size, split_number,base_path+"/splitted_files/" , checks)
    splitten_end = time.time()
    print "Time for splitting: " + str(time_difference(splitten_end - splitten_start))
    

    # 2. mappen
    # Listen der splitted files
    print "mapping... "
    
    array_job_ID, path_to_mapped_files, mapping_process = write_array_file(bwa_path, "/vol/biotools/bin/samtools", reference, forward_files[i], base_path+"/splitted_files/", mapping_threads, len(splitfilesList)/2)
    mappen_end = time.time()


    # 3. mergen
    
    mappedFilesList =  []
    for j in range(0,len(splitfilesList),2):
        print path_to_mapped_files
        print splitfilesList[j]
        print splitfilesList[j].replace('.fastq','.bam')
	mappedFilesList.append(os.path.join(path_to_mapped_files , splitfilesList[j].replace('.fastq','.bam')))

    print 'files mapped to:'
    print mappedFilesList
	    
    merged_bam = mergeBams("/vol/biotools/bin/samtools", forward_files[i] , mappedFilesList ,array_job_ID, mapping_process)
    print "merged bams to : " + merged_bam
    mergen_end = time.time()


    # 4. cleanup
    cleanup = False
    def purge(dir, pattern):
        for f in os.listdir(dir):
            if re.search(pattern, f):
                os.remove(os.path.join(dir, f))

    os.chdir(io_path)
    if testing:
        purge(io_path, os.path.split(final_file))
        if os.path.exists(base_path):
              print base_path + " exists"
              shutil.rmtree(base_path)
              print "removed"
              purge(error_log, "bwa_mapping*")
              purge(error_log, "bwa_merge*")
              purge(out_log, "bwa_mapping*")
              purge(out_log, "bwa_merge*")
    else:
                if os.path.exists(base_path+"/splitted_files") and cleanup == True:
                    print base_path+"/splitted_files exists"
                    shutil.rmtree(base_path+"/splitted_files")
                    print "removed"
                if os.path.exists(base_path+"/sam_files"):
                    print base_path+"/sam_files exists"
                    shutil.rmtree(base_path+"/sam_files")
                    print "removed"
                if os.path.exists(base_path+"/mapped_files"):
                    print base_path+"/mapped_files exists"
                    shutil.rmtree(base_path+"/mapped_files")
                    print "removed"
                if os.path.exists(base_path+"/sorted_mapped_files"):
                    print base_path+"/sorted_mapped_files exists"
                    shutil.rmtree(base_path+"/sorted_mapped_files")
                    print "removed"


    cleanup_end = time.time()


    dest = open(os.path.join(io_path, os.path.split(forward_files[i])[1].replace('_001.paired.fastq', "_" + str(split_size) + "s_" + str(split_number) + "n_" + str(mapping_threads) + "t.paired.txt")), 'w')

    stats = "Stats to " + forward_files[i] + ":\n"
    ss = "Split size: " + str(split_size) + "\n"
    sn = "Number of splits: " + str(split_number) + "\n"
    mt = "Mapping threads: " + str(mapping_threads) + "\n"
    ch = "Checks: " + str(checks) + "\n"
    init = "Init duration: " + time_difference(init_end - init_start) + "\n"
    splitten = "Splitting duration: " + time_difference(splitten_end - splitten_start) + "\n"
    mappen = "Mapping duration: " + time_difference(mappen_end - splitten_end) + "\n"
    mergen = "Merging duration: " + time_difference(mergen_end - mappen_end) + "\n"
    cleanup = "Clean up duration: " + time_difference(cleanup_end - mergen_end) + "\n"
    total = "Total duration: " + time_difference(cleanup_end - init_start) + "\n"

    dest.write(stats)
    dest.write(ss)
    dest.write(sn)
    dest.write(mt)
    dest.write(ch)
    dest.write(init)
    dest.write(splitten)
    dest.write(mappen)
    dest.write(mergen)
    dest.write(cleanup)
    dest.write(total)

    dest.close()


                    

    print "Job finished"

    


