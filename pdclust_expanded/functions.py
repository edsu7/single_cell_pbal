import functools
import os
import multiprocessing as mp
import glob
import time
import itertools
import sys
import pickle
import numpy as np
import scipy as sci
import pandas as pd
import subprocess
from pybedtools import BedTool
from plotly import tools

if len(sys.argv)<4:
	print("Warning not enough variables specified. Please supply a <OUT DIRECTORY> <INDEX FILE> <PROJECT NAME>")
	sys.exit()
elif len(sys.argv)==4:
	print ("number of jobs and threads not supplied. Defaulting to 4")
	threads=4
	jobs=4
elif len(sys.argv)<5:
	print ("number of jobs not supplied. Defaulting to 4")
	jobs=4
else:
	threads=int(sys.argv[4])
	jobs=int(sys.argv[5])

out_dir=sys.argv[1]
index_file=sys.argv[2]
project_name=sys.argv[3]
print("output directory:"+out_dir)
print("output index file:"+index_file)
#######################################
def calcFractionalMethylation(out_dir):
    """
    Function for converting out_dir/extract/**/*_cpg.bed.gz into strand collapsed out_dir/extract/**/*.fractional_methylation.bed.gz 
    """
    print("".join(["#"]*18))
    print("Generating fractional methylation calls")
    t0=time.time()
    reference_cpgs="/ref/CG.bed.gz"
    
    
    
    for x in glob.iglob(out_dir+"/extract/**/*_cpg.bed.gz", recursive=True):
        print("fractional Methylation "+x)
        fractionalMethylation(reference_cpgs,out_dir,x)
    
    print("Run time:"+str(time.time()-t0))

#######################################
def fractionalMethylation(reference_cpgs,out_dir,cpg_file):
    known_CpG=BedTool(reference_cpgs)
    known_CpG.map(BedTool.from_dataframe(pd.read_csv(cpg_file,compression='gzip',skiprows=0,
                        names=[
                        "chr",
                        "start",
                        "stop",
                        "name",
                        "score",
                        "strand",
                        "display_start",
                        "display_end",
                        "color",
                        "coverage",
                        "methylation",
                        "ref_geno",
                        "sample_geno",
                        "quality_score"
                        ],
                        usecols=[
                            "chr",
                            "start",
                            "stop",
                            "strand",
                            "coverage",
                            "methylation"
                        ],
                        sep='\t'
                       )\
            .query("coverage>0")\
            .assign(methylated = lambda row : round(row["coverage"]*row["methylation"]/100))\
            .assign(unmethylated = lambda row : row["coverage"]-row['methylated'])\
            .assign(start = lambda row : row["start"].astype(int))\
            .assign(stop = lambda row : row["stop"].astype(int))
        ),
        c=[7,8])\
    .to_dataframe()\
    .rename(columns={"name":"methylated","score":"unmethylated"})\
    .query("methylated!='.' and unmethylated!='.'")\
    .replace(".",0)\
    .assign(coverage = lambda row : row['methylated'].astype(int)+row['unmethylated'].astype(int))\
    .assign(frac_meth = lambda row : round(row['methylated'].astype(int)/row['coverage'],2))\
    .to_csv(cpg_file.replace("_cpg.bed.gz",".fractional_methylation.bed.gz"),compression='gzip',sep='\t',index=False,header=False)
    BedTool.cleanup()
    #.query("methylated!='.' or unmethylated!='.'")\
    #.query("coverage>=3")\
            
######################################
def distributeJobs(jobs,total_cmd_list):
    pool = mp.Pool(jobs,maxtasksperchild=1)  
    pool.map(runCommand,total_cmd_list)
    pool.close()

##########################################
def runCommand(sample_cmd_list):
    
    for cmds in sample_cmd_list:
        out_dir=cmds[0]
        cmd=cmds[1]
        cmd_name=cmds[2]
        variable=cmds[3]
        
        print(" ".join(cmd))

        subprocess.run(["mkdir","-p",out_dir+"/log"])
        if variable=='shell':
            f=open(cmd[cmd.index(">")+1:][0], "w")
            result=subprocess.run(cmd[:cmd.index(">")],capture_output=True)
            f.write(result.stdout.decode('utf-8'));f.close()        
        elif variable=='log':
            f=open(out_dir+"/log/"+cmd_name+".stdout", "w")
            result=subprocess.run(cmd,capture_output=True)
            f.write(result.stdout.decode('utf-8'));f.close()

            f=open(out_dir+"/log/"+cmd_name+".stderr", "w")
            f.write(result.stderr.decode('utf-8'));f.close()
        elif variable=='run':
            pass
            subprocess.run(cmd)
        else:
            pass
            
            
#######################################
def runGEMbs(out_dir):
    print("".join(["#"]*18))
    t0=time.time()
    cmd=["gemBS","prepare","-c",out_dir+"config.txt","-t",out_dir+"metadata.csv"]
    runCommand([[out_dir,cmd,"gemBS_prepare",'log']])
    
    cmd=["gemBS","index"]
    runCommand([[out_dir,cmd,"gemBS_index",'log']])

    
    cmd=["gemBS","--loglevel","debug","map"]
    runCommand([[out_dir,cmd,"gemBS_map",'log']])
    

    cmd=["gemBS","--loglevel","debug","merge-bams"]
    runCommand([[out_dir,cmd,"gemBS_merge",'log']])
    
    total_cmd_list=[]
    for x in glob.iglob(out_dir+"/mapping/**/*.bam", recursive=True):
        sample_cmd_list=[]
        cmd=[
        "java",
        "-Xmx4g",
        "-jar",
        '/usr/local/anaconda/share/picard-2.22.3-0/picard.jar',
        'MarkDuplicates',
        "I="+x,
        "O="+x.replace(".bam",".dups.marked.sorted.bam"),
        "M="+x.replace(".bam","_metrics.txt"),
        "VALIDATION_STRINGENCY=SILENT",
        "ASSUME_SORTED=true",
        "TMP_DIR="+out_dir+"/mapping"
        ]
        #runCommand(out_dir,cmd,"gemBS_index",'log')
        sample_cmd_list.append([out_dir,cmd,"markDup_"+x.split("/")[-1],'run'])
        
        cmd=["mv",
             x.replace(".bam",".dups.marked.sorted.bam"),
             x]
        sample_cmd_list.append([out_dir,cmd,"mv_"+x.split("/")[-1],'log'])
        #runCommand(out_dir,cmd,"markDup_"+x.split("/")[-1],'log')

        cmd=["md5sum",x,">",x.replace(".bam",".bam.md5")]
        sample_cmd_list.append([out_dir,cmd,"recalc_md5sum",'shell'])
        #runCommand(out_dir,cmd,"recalc_md5sum",'shell')

        cmd=["samtools",
             "flagstat",
             "-@"+str(threads),
             x.replace(".dups.marked.sorted.bam",".bam"),
             ">",
             x.replace(".dups.marked.sorted.bam",".bam").replace(".bam",".flagstat")]
        sample_cmd_list.append([out_dir,cmd,"flagstat",'shell'])
        #runCommand(out_dir,cmd,"flagstat",'shell')

        cmd=["samtools","index","-c","-@"+str(threads),x]
        sample_cmd_list.append([out_dir,cmd,"csi",'run'])
        #runCommand(out_dir,cmd,"csi",'run')
    
    distributeJobs(jobs,total_cmd_list)
        
    cmd=["gemBS","--loglevel","debug","call"]
    runCommand([[out_dir,cmd,"gemBS_call",'log']])

    cmd=["gemBS","--loglevel","debug","extract"]
    runCommand([[out_dir,cmd,"gemBS_extract",'log']])
    print("Run time:"+str(time.time()-t0))
    

#######################################
def runTrimGalore(csv_file,out_dir):
    print("".join(["#"]*18))
    t0=time.time()
    tmp=pd.read_csv(csv_file,sep='\t',names=['index','read1','read2'])
    total_cmd_list=[]
    for x in tmp.values.tolist():
        cmd=["trim_galore","--clip_R1","6","--clip_R2","6","--paired",x[1],x[2],"-o",out_dir+"/fastq/","-j",str(threads),"--gzip","--fastqc"]
        total_cmd_list.append([[out_dir,cmd,"trim_"+x[0],"log"]])
    
    distributeJobs(jobs,total_cmd_list)
    
    print("Run time:"+str(time.time()-t0))

#######################################
def setUpMetadata(csv_file,project_name,out_dir):
    print("".join(["#"]*18))
    print("SETTING UP metadata.csv")
    tmp=pd.read_csv(csv_file,sep='\t',names=['index','read1','read2'])\
    .assign(end_1 = lambda row : out_dir+"/fastq/"+row['read1'].str.split("/").str[-1].str.replace(".fastq.gz","_val_1.fq.gz"))\
    .assign(end_2 = lambda row : out_dir+"/fastq/"+row['read2'].str.split("/").str[-1].str.replace(".fastq.gz","_val_2.fq.gz"))\
    .assign(Barcode = lambda row : row['index'])\
    .assign(Library = lambda row : project_name)\
    .assign(file_id = lambda row : row['Library']+"_"+row['index'])\
    .assign(file_name = lambda row : row['Library']+"_"+row['index'])\
    .loc[:,["Barcode","Library","file_id","end_1","end_2","file_name"]]
    tmp.to_csv(out_dir+"/metadata.csv",sep=',',index=False)
    return(out_dir+"/metadata.csv")
########################################
def readSingleCellIndex(index_file_path):
    print("".join(["#"]*18))
    csv_file=None
    read_error=False
    if os.path.exists(index_file_path):
        try:
            csv_file=pd.read_csv(index_file_path,sep='\t',names=["index","read1","read2"])
            if (csv_file.shape[0]>1) and (csv_file.shape[1]>=3):
                print("Format looks good. Continuing")
            else:
                print("Format Error")
                sys.exit(1)
        except:
            print("Problem openning {}. Check file path or if file is in appropriate format".format(index_file_path))
    else:
        print("{} does not exist".format(index_file_path))
        sys.exit(1)  
    
    for x in csv_file.read1.values.tolist()+csv_file.read2.values.tolist():
        if not os.path.exists(x):
            print("{} does not exist".format(x))
            read_error=True
    
    if read_error:
        print("Issues persist. Please address them and run again")
        sys.exit(1)
        
    return(csv_file)


############################################
def gemBS_ConfigurationSetup(output,project_name):
    print("".join(["#"]*18))
    print("SETITNG UP config.txt")
    f=open(output+"/config.txt","w+")
    ### RUNNING PARAMETERS
    ### AS PER IHEC STANDARDS ; SEE GITHUB IF CHANGES ARE NECESSARY
    ### WORKING DIRECTORY - mounted ###
    f.write(
    "base="+output+" ### if mounted by following example do not change ###"+"\n"+\
    ""+"\n"+\
    "sequence_dir = ${base}/fastq/@SAMPLE    # @SAMPLE and @BARCODE are special"+"\n"+\
    "bam_dir = ${base}/mapping/@BARCODE      # variables that are replaced with"+"\n"+\
    "bcf_dir = ${base}/calls/@BARCODE        # the sample name or barcode being"+"\n"+\
    "extract_dir = ${base}/extract/@BARCODE  # worked on during gemBS operation"+"\n"+\
    "report_dir = ${base}/report"+"\n"+\
    ""+"\n"+\
    "### REFERENCES - mounted ###"+"\n"+\
    "### if mounted by following example do not change ###"+"\n"+\
    "reference = /ref/hg38_no_alt.fa"+"\n"+\
    "index_dir = /ref"+"\n"+\
    "extra_references = /ref/conversion_control.fa"+"\n"+\
    ""+"\n"+\
    "# General project info"+"\n"+\
    "project = "+project_name+" ### SPECIFIC PROJECT TITLE ###"+"\n"+\
    "species = hg38"+"\n"+\
    ""+"\n"+\
    "# Default parameters"+"\n"+\
    "threads = "+str(threads)+"\n"+\
    "jobs = "+str(jobs)+" ### MODIFY FOR NEED BE ###"+"\n"+\
    ""+"\n"+\
    "[index]"+"\n"+\
    ""+"\n"+\
    "sampling_rate = 4"+"\n"+\
    ""+"\n"+\
    "[mapping]"+"\n"+\
    ""+"\n"+\
    "non_stranded = True ### TOGGLE TO TRUE FOR PBAL ###"+"\n"+\
    "remove_individual_bams = True"+"\n"+\
    "underconversion_sequence = NC_001416.1"+"\n"+\
    "overconversion_sequence = V01146.1"+"\n"+\
    ""+"\n"+\
    "[calling]"+"\n"+\
    ""+"\n"+\
    "mapq_threshold = 10"+"\n"+\
    "qual_threshold = 13"+"\n"+\
    "reference_bias = 2"+"\n"+\
    "left_trim = 0"+"\n"+\
    "right_trim = 0"+"\n"+\
    "keep_improper_pairs = True ### TOGGLE TO TRUE FOR PBAL ###"+"\n"+\
    "keep_duplicates = False ### TOGGLE TO TRUE FOR RRBS -  ###"+"\n"+\
    "haploid = False"+"\n"+\
    "conversion = auto"+"\n"+\
    "remove_individual_bcfs = True"+"\n"+\
    "contig_pool_limit = 25000000"+"\n"+\
    ""+"\n"+\
    "[extract] # extract specific section"+"\n"+\
    ""+"\n"+\
    "strand_specific = True"+"\n"+\
    "phred_threshold = 10"+"\n"+\
    "make_cpg = True"+"\n"+\
    "make_non_cpg = False"+"\n"+\
    "make_bedmethyl = True"+"\n"+\
    "make_bigwig = True"+"\n"
    )
    f.close()
############################################
def freec_ConfigurationSetup(output,project_name):
    print("".join(["#"]*18))
    print("SETITNG UP freec configs")            
    cmd=['mkdir','-p',out_dir+"/cnv"]
    runCommand([[out_dir,cmd,"making_cnv_dir","run"]])
            
    for x in glob.iglob(out_dir+"/mapping/**/*.bam", recursive=True):
        name=x.split("/")[-1].replace(".bam","")
        cmd=['mkdir','-p',out_dir+"/cnv/"+name]
        runCommand([[out_dir,cmd,"mkdir_cnv","run"]])
        
        f=open(out_dir+"/cnv/"+name+"/config.txt","w+")
        f.write(
        "[general]"+"\n"+\
        "chrFiles=/ref2"+"\n"+\
        "chrLenFile=/ref/hg38_no_alt_autosomes.contig.sizes"+"\n"+\
        "maxThreads="+str(threads)+"\n"+\
        "ploidy=2"+"\n"+\
        "samtools=//usr/local/anaconda/bin/samtools"+"\n"+\
        "window=5000000"+"\n"+\
        "telocentromeric=5000000"+"\n"+\
        "outputDir="+out_dir+"/cnv/"+name+"\n"+\
        "sex=XY"+"\n"+\
        "minExpectedGC=0.39"+"\n"+\
        "maxExpectedGC=0.51"+"\n"+\
        "\n"+\
        "[sample]"+"\n"+\
        "mateFile="+out_dir+"/cnv/"+name+"/"+name+".dedup.bam"+"\n"+\
        "inputFormat=BAM\n"
        )
        f.close()
    
############################################
def runControlFREEC(out_dir):
    """
    Function for copy number caller
    """
    print("".join(["#"]*18))
    t0=time.time()   
    total_cmd_list=[]
    for x in glob.iglob(out_dir+"/mapping/**/*.bam", recursive=True):
        sample_cmd_list=[]
        name=x.split("/")[-1].replace(".bam","")
        cmd=["samtools",
             "view",
             x,
             "-@"+str(threads),
             "-h",
             "-b",
             "-F516",
             "-o",
             x.replace(".bam",".dedup.bam").replace("mapping","cnv")
        ]

        sample_cmd_list.append([out_dir,cmd,"dedup_"+name,"run"])
        cmd=["freec","-conf",out_dir+"/cnv/"+name+"/config.txt"]

        sample_cmd_list.append([out_dir,cmd,"freec_"+name,"log"])

        cmd=["rm",
             x.replace(".bam",".dedup.bam").replace("mapping","cnv")
        ]
        sample_cmd_list.append([out_dir,cmd,"rm_"+name+"_dedup","run"])
        total_cmd_list.append(sample_cmd_list)
    
    distributeJobs(jobs,total_cmd_list)
    print("Run time:"+str(time.time()-t0))
############################################
        

readSingleCellIndex(index_file)
gemBS_ConfigurationSetup(out_dir,project_name)
setUpMetadata(index_file,project_name,out_dir)
runTrimGalore(index_file,out_dir)
runGEMbs(out_dir)
calcFractionalMethylation(out_dir)
freec_ConfigurationSetup(out_dir,project_name)
runControlFREEC(out_dir)
print("".join(["#"]*18)+"\nFinished")
