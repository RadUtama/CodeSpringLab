
import config
import pandas as pd
import os
import time
import re
from IPython.display import IFrame,clear_output
import shutil
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame

project_name=config.project_name

def Tree():
    
    cmd = "tree ../../csl_results/"
    print(os.popen(cmd).read())

def DeleteJobs(jobid):
    
    for i in range(len(jobid)):
        command = "qdel "+jobid[i]
        jobdel=os.popen(command).read().splitlines()
        print('Deleting job ID: '+jobid[i])
        while True:
            command = "qstat"
            del_status=os.popen(command).read().splitlines()
            if del_status == []:
                break
            else:
                del_stat = pd.DataFrame(del_status)[0][2:].str.split(expand=True)[0].str.contains(jobid[i]).any()
                if del_stat == False:
                    break
            time.sleep(5)
            
def DeleteOneJob(jobid):
    
    command = "qdel "+jobid[0]
    jobdel=os.popen(command).read().splitlines()
    print('Deleting job ID: '+jobid[0])
    while True:
        command = "qstat"
        del_status=os.popen(command).read().splitlines()
        if del_status == []:
            break
        else:
            del_stat = pd.DataFrame(del_status)[0][2:].str.split(expand=True)[0].str.contains(jobid[0]).any()
            if del_stat == False:
                break
        time.sleep(5)
            
def Qstat(jobid):
        
    print("==================================")
    starttime = time.time()
    while True:
        command = "qstat"
        running_status = os.popen(command).read().splitlines()
    
        clear_output(wait=True)
    
        if running_status == []:
            print("All jobs done !" )
            print(".........................")
            print("Running time: "+str( round( ((time.time() - starttime) / 60) , 2 ) )+" minutes")
            runstat = "done"
        elif pd.DataFrame(running_status)[0][2:].str.split(expand=True)[0].str.contains('|'.join(jobid)).any():
            
            all_stat = pd.DataFrame(running_status)[0][2:].str.split(expand=True)[4]
            id_sub = pd.DataFrame(running_status)[0][2:].str.split(expand=True)[0].str.contains('|'.join(jobid))
            id_stat = all_stat[id_sub]
            status_stat = id_stat.value_counts().rename(index={'r':'Number of jobs still running = ',
                                                    'qw':'Number of jobs waiting in line = ',
                                                    't':'Number of jobs about to run = ',
                                                    'Rt':'Number of jobs about to run = ',
                                                    'Rr':'Number of jobs about to run = ',
                                                    'Eqw':'Number of jobs cannot run (server error) = ',
                                                   })
            print(status_stat)
            print(".........................")
            print("Running time: "+str( round( ((time.time() - starttime) / 60) , 2 ) )+" minutes")
            runstat = "undone"
        else:
            print("All jobs done !" )
            print(".........................")
            print("Running time: "+str( round( ((time.time() - starttime) / 60) , 2 ) )+" minutes")
            runstat = "done"
    
        if runstat=="done":
            break
        time.sleep(10)

def ListDir(directory):
        
    print("Here's the list of contents:")
    print("Index")
    dirlist = pd.Series(os.listdir(directory))
    print(dirlist)
    
    return dirlist

def filetransfer_Prep():
        
    global project_name
    os.makedirs("../../csl_results/"+project_name+"/data/",exist_ok=True)
    #os.makedirs("../../csl_results/"+project_name+"/data/manifest/",exist_ok=True)
    os.makedirs("../../csl_results/"+project_name+"/log/",exist_ok=True)
    
    if os.path.exists("../../csl_results/"+project_name+"/log/output_listFastq.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_listFastq.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_listFastq.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_listFastq.txt")
   
    read_path_destination = "../../csl_results/"+project_name+"/data/fastq/"
    print("Read files will be copied to ../../csl_results/"+project_name+"/data/fastq/")
    print("==================================")
    print("Here's the list of available genomes:")
    genome_list = pd.Series(['human','mouse'])
    print("Index")
    print(genome_list)
    print("==================================")
    print("Specify the index to the genome:(e.g 0)")
    print("\033[91m"+"If you want to use our example dataset, type the number"+"\033[94m"+" 0"+"\x1b[0m")
    genome_index = int(input())
    genome = genome_list[genome_index]
    print("==================================")
    print("Copy the path to your original read files folder:")
    print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
    print("../scripts_DoNotTouch/test/fastq_atac/")
    read_path_original = input()
    print("==================================")
    print("Copy the path to design matrix folder (If it's in your home folder, type tilde sign ~):")
    print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
    print("../scripts_DoNotTouch/test/manifest_atac/")
    inpath_design = input()
    print("==================================")
    print("You'll be working with "+"\033[91m"+project_name+"\x1b[0m"+" folder")
    print("Re-running any cell will overwrite exisiting outputs in folder "+"\033[91m"+project_name+"\x1b[0m")
    print("If you don't want to overwrite, please re-run this cell and specify different unique"+"\033[91m project_name")
    
    scriptpath_listdir = "../scripts_DoNotTouch/fastq/qsub_listdir.sh"
    scriptpath_copy = "../scripts_DoNotTouch/fastq/qsub_copy.sh"
    
    des=pd.read_table(inpath_design+"/design_matrix.txt")
    filename=des.iloc[:,len(des.columns)-1]
    if filename.str.contains('_R2_').any() ==True:
        pairing = 'y'
    else:
        pairing = 'n'

    #command = "source "+scriptpath_listdir+" "+read_path_original+" "+project_name
    #joblist=os.popen(command).read().splitlines()
    #print(joblist)
    
    #return read_path_original,read_path_destination,scriptpath_copy,scriptpath_listdir,genome,pairing,inpath_design
    return read_path_original+"/",read_path_destination+"/",scriptpath_copy,genome,pairing,inpath_design+"/"

def filetransfer_ListDir(read_path_original):
    
    global project_name
    dirfileset = ['empty']
    
    if os.path.exists("../../csl_results/"+project_name+"/log/output_copyFastq.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_copyFastq.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_copyFastq.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_copyFastq.txt")
    
    try:
        print("Here's the list of files in the original folder:")
        print("Index")
        listori = pd.read_csv("../../csl_results/"+project_name+"/log/output_listFastq.txt",header=None)
        listori = listori[listori[0].str.endswith('fastq.gz')]
        print(listori)

        dirfileset = read_path_original + listori[listori[0].str.endswith('fastq.gz')]
    except OSError:
        print("Access to view files in original directory is still pending")
    
    return dirfileset[0]

def filetransfer_PrepDirect():
    
    print("========================================")
    print("Specify the path to fastq folder used for QC:")
    read_path_destination = input()
    print("========================================")
    
    return read_path_destination+"/"

def filetransfer_Copy(read_path_original,scriptpath_copy):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/data/fastq") :
        shutil.rmtree("../../csl_results/"+project_name+"/data/fastq")
    
    jobid = []

    command = "source "+scriptpath_copy+" "+read_path_original+" "+"../../csl_results/"+project_name+"/data/fastq"+" "+project_name
    job = os.popen(command).read().splitlines()
    print(job)
    jobid.append(job[0].split(' ')[2])
    
    return jobid

def filetransfer_ListDest(directory):
    
    global project_name
    
    print("Here's the list of contents:")
    print("Index")
    dirlist = pd.Series(os.listdir(directory))
    
    dirlist = dirlist[dirlist.str.endswith('fastq.gz')]
    print(dirlist)
    
    dirfileset = directory + dirlist
    
    rmhidden = [shutil.rmtree(f) for f in os.listdir("../../csl_results/"+project_name+"/data/fastq") if f.startswith(".")]
    
    return dirfileset

def filetransfer_Convert(directory,inpath_design):
    
    des=pd.read_table(inpath_design+"/design_matrix.txt")
    filename=des.iloc[:,len(des.columns)-1]

    dirlist = pd.Series(os.listdir(directory))    
    dirlist = dirlist[dirlist.str.endswith('fastq.gz')]
    for i in range(len(dirlist)):
        for j in range(len(filename)):
            if dirlist.iloc[i] in filename.iloc[j]:
                newname=des.iloc[j,0]+re.sub(r'^.*?_R', '_R', dirlist.iloc[i])
                os.rename(directory+dirlist.iloc[i],directory+newname)
     
    print("Here's the list of name-converted read files:")
    print("Index")
    dirlist = pd.Series(os.listdir(directory))
    
    dirlist = dirlist[dirlist.str.endswith('fastq.gz')]
    print(dirlist)
    
    dirfileset = directory + dirlist
    
    return dirfileset

def fastqc_Prep(directory):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_fastQC.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_fastQC.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_fastQC.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_fastQC.txt")
    
    folder_fastqc = "fastqc"
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        directory = "../../csl_results/"+project_name+"/data/cutadapt/"
        folder_fastqc = "fastqc_cutadapt"
    print("========================================")
    
    readlist = pd.Series(os.listdir(directory))
    readlist = readlist[readlist.str.endswith('fastq.gz')]
    
    outdir_fastqc = "../../csl_results/"+project_name+"/data/"+folder_fastqc+"/"
    os.makedirs(outdir_fastqc,exist_ok=True)
    
    print("FastQC results will be stored in ../../csl_results/"+project_name+"/data/"+folder_fastqc+"/")
    
    scriptpath_fastqc = '../scripts_DoNotTouch/FastQC/qsub_fastqc.sh'

    return readlist,directory,outdir_fastqc,scriptpath_fastqc

def fastqc_PrepDirect():
    
    print("========================================")
    print("Specify the path to fastq folder used for QC:")
    read_path_destination = input()
    print("========================================")
    
    return read_path_destination+"/"

def fastqc_RunQC(readlist,outdir_fastqc,read_path_destination,scriptpath_fastqc):
            
    global project_name
    
    jobid = []
    for file in readlist:
        command = "source "+scriptpath_fastqc+" "+read_path_destination+file+" "+outdir_fastqc+"/."+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job)
        
        jobid.append(job[0].split(' ')[2])
    
    return jobid

def fastqc_ListDir(outdir_fastqc):
    
    dirlist = pd.Series(os.listdir(outdir_fastqc))
    dirlist = dirlist[dirlist.str.endswith('.html')]
    dirlist.index = range(len(dirlist))
    
    DataFrame(dirlist).to_html(outdir_fastqc+'/QC_list.html')
    dir_html = IFrame(outdir_fastqc+'/QC_list.html', width=1000, height=800)
    
    return dir_html

def fastqc_Visualization(outdir_fastqc):

    dirlist = pd.Series(os.listdir(outdir_fastqc))
    dirlist = dirlist[dirlist.str.endswith('.html')]
    dirlist.index = range(len(dirlist))
    #print(dirlist)
    
    print("Specify index to visualize HTML file:(e.g 0)")
    index_files = int(input())    
    
    file = dirlist[index_files]
    qc = IFrame(outdir_fastqc+file, width=1000, height=800)
    
    return qc

def cutadapt_Prep(directory,pairing):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_Cutadapt.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_Cutadapt.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_Cutadapt.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_Cutadapt.txt")
    
    readlist = pd.Series(os.listdir(directory))
    readlist = readlist[readlist.str.endswith('fastq.gz')]
    
    prefix = readlist.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    outdir_cutadapt = "../../csl_results/"+project_name+"/data/cutadapt/"
    os.makedirs(outdir_cutadapt,exist_ok=True)
    
    print("Trimmed reads results will be stored in ../../csl_results/"+project_name+"/data/cutadapt/")
    
    if pairing == "y":
        scriptpath_cutadapt = '../scripts_DoNotTouch/cutadapt_PE/qsub_cutadapt_PE.sh'
        read1_list = directory+'/'+prefix+'_R1_001.fastq.gz'
        read2_list = directory +'/'+prefix+'_R2_001.fastq.gz'
        trimmed1_list = outdir_cutadapt+'/'+prefix+'_R1_001.fastq.gz'
        trimmed2_list = outdir_cutadapt+'/'+prefix+'_R2_001.fastq.gz'
    else:
        scriptpath_cutadapt = '../scripts_DoNotTouch/cutadapt_SE/qsub_cutadapt_SE.sh'
        read1_list = directory+'/'+prefix+'_R1_001.fastq.gz'
        read2_list = directory +'/'+prefix+'_R2_001.fastq.gz'
        trimmed1_list = outdir_cutadapt+'/'+prefix+'_R1_001.fastq.gz'
        trimmed2_list = outdir_cutadapt+'/'+prefix+'_R2_001.fastq.gz'
        
    print("========================================")
    print("Specify the adapter for R1/read1:")
    print("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA for Illumina Universal TruSeq RNA")
    print("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC for Nextera Transposase ATAC")
    print("TGGAATTCTCGG for Illumina Small RNA 3' ")
    print("GATCGTCGGACT for Illumina Small RNA 5' ")
    adapter = input()
    print("========================================")
    print("Specify the adapter for R2/read2:")
    print("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT for Illumina Universal TruSeq RNA")
    print("CTGTCTCTTATACACATCTGACGCTGCCGACGA for Nextera Transposase ATAC")
    print("TGGAATTCTCGG for Illumina Small RNA 3' ")
    print("GATCGTCGGACT for Illumina Small RNA 5' ")
    adapter2 = input()
    print("========================================")
    print("Specify minimum length after trimming (default 20):")
    minlen = input()

    return adapter,adapter2,minlen,read1_list,read2_list,trimmed1_list,trimmed2_list,outdir_cutadapt,scriptpath_cutadapt

def cutadapt_PrepDirect():
    
    print("========================================")
    print("Specify the path to fastq folder used for adapter trimming:")
    read_path_destination = input()
    print("========================================")
    
    return read_path_destination+"/"

def cutadapt_RunTrimming(adapter,adapter2,minlen,read1_list,read2_list,trimmed1_list,trimmed2_list,scriptpath_cutadapt):
            
    global project_name
    
    jobid = []
    for i in range(len(read1_list)):
        command = "source "+scriptpath_cutadapt+" "+minlen+" "+adapter+" "+adapter2+" "+trimmed1_list[i]+" "+trimmed2_list[i]+" "+read1_list[i]+" "+read2_list[i]+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job)
        
        jobid.append(job[0].split(' ')[2])
    
    return jobid

def bowtie2_Prep(genome,pairing,read_dir):
        
    global project_name
    
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        read_dir = "../../csl_results/"+project_name+"/data/cutadapt/"
    print("========================================")
    
    prefix = pd.Series(os.listdir(read_dir))
    prefix = prefix[prefix.str.endswith('fastq.gz')]
    prefix = prefix.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    out_dir = "../../csl_results/"+project_name+"/data/bowtie2/"
   
    for i in range(len(prefix)):
        os.makedirs(out_dir+prefix[i],exist_ok=True)

    print("Bowtie2 alignment results will be stored in ../../csl_results/"+project_name+"/data/bowtie2/")
    
    if genome == 'mouse':
        genome_index_path = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/GRCm39_M29_gencode_bowtie2index/GRCm39_M29_gencode"
    elif genome == 'human':
        genome_index_path = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/hg38_p13_gencode_bowtie2index/hg38_p13_gencode"
    
    read1_list = read_dir+'/'+prefix+'_R1_001.fastq.gz'
    read2_list = read_dir+'/'+prefix+'_R2_001.fastq.gz'
    out_prefix_list = out_dir+prefix+'/'+prefix
    
    if pairing == 'y':
        scriptpath_bowtie2 = '../scripts_DoNotTouch/bowtie2/qsub_bowtie2_PE.sh'
    else:
        scriptpath_bowtie2 = '../scripts_DoNotTouch/bowtie2/qsub_bowtie2_SE.sh'

    return genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,scriptpath_bowtie2

def bowtie2_PrepDirect():
    
    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Are the reads paired-end:(e.g y/n)")
    pairing = input()
    print("========================================")
    print("Specify the path to fastq folder used for alignment:")
    read_path_destination = input()
    print("========================================")
    
    return genome,pairing,read_path_destination+"/"

def bowtie2_RunAlignment(genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,scriptpath_bowtie2):
        
    global project_name
    
    cmd_dir = "cat > ../../csl_results/"+project_name+"/log/error_bowtie2.txt"
    cmd_dir_run = os.popen(cmd_dir)
    
    jobid = []
    for i in range(len(out_prefix_list)):
        command = "source "+scriptpath_bowtie2+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        job = os.popen(command).read().strip().splitlines()
        #job = os.popen(command).read().splitlines()
        print(job)
        jobid.append(job[0].split(' ')[2])
    
    return jobid

def bowtie2_ListDir(directory):
    
    global project_name
    bowtie2logdir = "../../csl_results/"+project_name+"/data/bowtie2_summary/"
    os.makedirs(bowtie2logdir,exist_ok=True)
    
    dirlist = pd.Series(os.listdir(directory))
    
    i = 0

    for file in dirlist:
        i += 1
        
        logfile = directory+file+"/"+file+"Log.final.out"

        log_df = pd.read_csv(logfile,skiprows=15,nrows=15)
        log_raw = log_df.rename({log_df.columns[0]:file},axis='columns')
        
        if i == 1 :
            log_matrix = log_raw
        else :
            log_matrix = pd.concat([log_matrix,log_raw],axis=1)

    log_matrix.to_csv(bowtie2logdir+'summary_matrix.txt',sep='\t')        
    
    return log_matrix

    
def bowtie2_Visualization(out_dir):

    dirlist = DataFrame(pd.Series(os.listdir(out_dir)))
    dirlist.index = range(len(dirlist))
    print(dirlist)
    
    print("========================================")
    print("Specify index to visualize insert size frequency:(e.g 0)")
    index_files = int(input())    
    
    file = dirlist.loc[index_files,0]
    qc = IFrame(out_dir+"/"+file+"/"+file+"_insert_size_histogram-1.jpg", width=1000, height=800)
    
    return qc
    

def macs2_ListDir(prefix,count_prefix_list):
    
    global project_name
    outpath_macs2 = "../../csl_results/"+project_name+"/data/macs2/"
    os.makedirs(outpath_macs2,exist_ok=True)
    print("MACS2 summary matrix is stored in ../../csl_results/"+project_name+"/data/macs2/")
    
    i = 0

    for file in count_prefix_list:
        i += 1
        
        logfile = file+"_counts.txt.summary"

        log_df = pd.read_table(logfile,comment='#',header=None,sep="\t",index_col=[0])
        log_raw = log_df.rename({log_df.columns[0]:prefix[i-1]},axis='columns')
        if i == 1 :
            log_matrix = log_raw
        else :
            log_matrix = pd.concat([log_matrix,log_raw],axis=1)

    log_matrix.to_csv(outpath_counts+'featurecounts_summary.txt',sep='\t')        
    
    return log_matrix

def macs2_Prep(genome,out_dir,pairing):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_macs2.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_macs2.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_macs2.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_macs2.txt")
    
    prefix = pd.Series(os.listdir(os.path.expanduser(out_dir)))
    
    macs2_dir = "../../csl_results/"+project_name+"/data/macs2/"
   
    for i in range(len(prefix)):
        os.makedirs(macs2_dir+prefix[i],exist_ok=True)
    
    if genome == 'mouse':
        genomesize = "1.87e+9"
        chromsize = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/chrom.sizes"
        anno_onlyChrNoMito = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_onlyChrNoMito.bed"
    elif genome == 'human':
        genomesize = "2.7e+9"
        chromsize = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/chrom.sizes"
        anno_onlyChrNoMito = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_onlyChrNoMito.bed"
    
    if pairing == "y":
        scriptpath_macs2 = '../scripts_DoNotTouch/MACS2/qsub_macs2_PE.sh'
    else:
        scriptpath_macs2 = '../scripts_DoNotTouch/MACS2/qsub_macs2_SE.sh'

    print("Remove all kinds of read duplicates:(e.g y/n)")
    removeDup = input()
    if removeDup == "y":
        bed_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByName_removeDup.out.bed'
    else:
        bed_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByName.out.bed'
    
    macs2_prefix_list = macs2_dir+prefix+'/'

    return scriptpath_macs2,genomesize,chromsize,bed_list,macs2_prefix_list,prefix,anno_onlyChrNoMito

def macs2_PrepDirect():
    
    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Are the reads paired-end:(e.g y/n)")
    pairing = input()
    print("========================================")
    print("Specify the path to alignment folder used for peak calling:")
    out_dir = input()
    print("========================================")
    
    return genome,pairing,out_dir+"/"

def macs2_RunPeakCalling(scriptpath_macs2,genomesize,chromsize,bed_list,macs2_prefix_list,prefix,anno_onlyChrNoMito):
     
    global project_name
    
    jobid = []
    for i in range(len(bed_list)):
        command = "source "+scriptpath_macs2+" "+prefix[i]+" "+bed_list[i]+" "+genomesize+" "+chromsize+" "+macs2_prefix_list[i]+" "+anno_onlyChrNoMito+" "+project_name
        job = os.popen(command).read().strip().splitlines()
        #job = os.popen(command).read().splitlines()
        print(job)
        jobid.append(job[0].split(' ')[2])
    
    return jobid

def bowtie2_Prep(inpath_design):
        
    global project_name
    
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        read_dir = "../../csl_results/"+project_name+"/data/cutadapt/"
    print("========================================")
    
    prefix = pd.Series(os.listdir(read_dir))
    prefix = prefix[prefix.str.endswith('fastq.gz')]
    prefix = prefix.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    out_dir = "../../csl_results/"+project_name+"/data/bowtie2/"
   
    for i in range(len(prefix)):
        os.makedirs(out_dir+prefix[i],exist_ok=True)

    print("Bowtie2 alignment results will be stored in ../../csl_results/"+project_name+"/data/bowtie2/")
    
    read1_list = read_dir+'/'+prefix+'_R1_001.fastq.gz'
    read2_list = read_dir+'/'+prefix+'_R2_001.fastq.gz'
    out_prefix_list = out_dir+prefix+'/'+prefix
    
    scriptpath_homer_tag = '../scripts_DoNotTouch/Homer/qsub_homer_tag.sh'

    return genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,scriptpath_bowtie2


def featurecounts_CreateCountMatrix():

    global project_name
    inpath_counts = "../../csl_results/"+project_name+"/data/featurecounts/"
    outpath_counts = "../../csl_results/"+project_name+"/data/counts/"
    print("Count matrix is stored in ../../csl_results/"+project_name+"/data/counts/")

    filelist=sorted([f for f in os.listdir(inpath_counts) if not f.startswith('.')])

    i = 0

    for file in filelist:
        i += 1
        count_df = pd.read_table(os.path.join(inpath_counts,file,file+'_counts.txt'),comment='#',header=[0],index_col=[0])
        count_raw = count_df.drop(['Chr','Start','End','Strand','Length'],axis=1)
        count_raw = count_raw.rename({count_raw.columns[0]:file},axis='columns')
        if i == 1 :
            count_matrix = count_raw
        else :
            count_matrix = pd.concat([count_matrix,count_raw],axis=1)

    count_matrix.columns=count_matrix.columns.str.rstrip('_counts.txt')
    count_matrix.to_csv(outpath_counts+'count_matrix.txt',sep='\t')

    return outpath_counts,count_matrix
    
def deseq2_Prep(inpath_design):
    
    global project_name
    outpath = "../../csl_results/"+project_name+"/data/deseq2/"
    os.makedirs(outpath,exist_ok=True)
    print("DESeq2 results are stored in ../../csl_results/"+project_name+"/data/deseq2/")

    design = pd.read_table(inpath_design+'/design_matrix.txt',index_col=0)
    design = design.iloc[:,:len(design.columns)-1]
    
    print("========================================")
    print("Here's the list of phenotypes/conditions/experiments")
    design_var=[]
    for i in range(len(design.columns)):
        design_var.append(design.columns[i])
        print(design.columns[i]+':')
        print(set(design.iloc[:,i]))
    
    print("========================================")
    print("Which phenotype/condition/replicate/batch should be the reference/baseline?(e.g control)")
    refcond = input()
    
    print("========================================")
    print("Which phenotype/condition/replicate/batch to compare?(e.g treated)")
    compared = input()
    
    scriptpath_deseq2 = '../scripts_DoNotTouch/DESeq2/qsub_deseq2.sh'
    Rpath_deseq2 = '../scripts_DoNotTouch/DESeq2/DESeq2.R'
    
    return scriptpath_deseq2,Rpath_deseq2,outpath,refcond,compared,design_var

def deseq2_PrepDirect():
    
    print("========================================")
    print("Specify the path to folder containing count_matrix.txt used for DE:")
    outpath_counts = input()
    print("========================================")
    print("Specify the path to folder containing design_matrix.txt used for DE:")
    inpath_design = input()
    print("========================================")
    
    return outpath_counts+"/",inpath_design+"/"

def deseq2_RunDE(scriptpath_deseq2,Rpath_deseq2,inpath_counts,inpath_design,outpath,refcond,compared):
    
    global project_name
    
    jobid = []
    command = "source "+scriptpath_deseq2+" "+Rpath_deseq2+" "+inpath_counts+"/count_matrix.txt"+" "+inpath_design+"/design_matrix.txt"+" "+outpath+" "+refcond+" "+compared+" "+project_name
    job = os.popen(command).read().splitlines()
    print(job)
    jobid.append(job[0].split(' ')[2])

    return jobid
    
def visualization_PrepDirect():
    
    print("========================================")
    print("Specify the path to folder containing design_matrix.txt used for DE:")
    inpath_design = input()
    print("========================================")
    print("Specify the path to folder containing DESeq2 results:")
    outpath = input()
    print("========================================")
    print("Which phenotype/condition/replicate/batch should be the reference/baseline?(e.g control)")
    refcond = input()
    print("========================================")
    print("Which phenotype/condition/replicate/batch to compare?(e.g treated)")
    compared = input()
    
    return inpath_design+"/",outpath+"/",refcond,compared
    
def visualization_heatmap(inpath_design,outpath,refcond,compared):

    print("Provide path to folder containing genelist.txt (max 50 genes). To plot the top 50 differential genes from DESeq2 instead, type 'top50'")
    genepath = input()
    if genepath == "top50":
        genelist = pd.read_table(outpath+'/DEG_'+compared+'_vs_'+refcond+'(ref).txt',index_col=0)
        genes = genelist.index[:50]
    else:
        genelist = pd.read_table(genepath+'/genelist.txt',index_col=0,header=None)
        genes = genelist.index[:50]
    
    design = pd.read_table(inpath_design+"/design_matrix.txt",index_col=0)
    design = design.iloc[:,:len(design.columns)-1]
    count_norm = pd.read_table(outpath+'/normalized_counts.txt',index_col=0)
    count_norm = count_norm.drop(['DESCRIPTION'],axis=1)
    count_norm_sig = count_norm[count_norm.index.isin(genes)]
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300 
    for i in range(len(design.columns)):
        lut = dict(zip(design.iloc[:,i].unique(), sns.color_palette("Pastel2")))
        col_colors = design.iloc[:,i].map(lut)
        heat = sns.clustermap(count_norm_sig,z_score=0,cmap='vlag',col_colors=col_colors)
        heat.savefig(outpath+"/heatmap_"+design.columns[i]+".png")

    return heat

def visualization_peakHeatmap(design_var):

    global project_name
    outpath = "../../csl_results/"+project_name+"/data/deseq2/"
    
    print("Here's the list of phenotype/condition/replicate/batch:")
    print("Index")
    design_var = pd.DataFrame(design_var)
    print(design_var)
    print("==================================")
    print("Which index of phenotype/condition/replicate/batch to view PC plot?")
    inpca = int(input())
    
    return outpath+"/",inpca,design_var
