
import config
import pandas as pd
import os
import time
import re
from IPython.display import IFrame,clear_output,HTML,Image
import shutil
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame
import gseapy as gp
from gseapy import GSEA,Biomart,dotplot,heatmap
import imgkit

project_name=config.project_name
param=config.parameters_exist

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
    global param
    os.makedirs("../../csl_results/"+project_name+"/data/",exist_ok=True)
    #os.makedirs("../../csl_results/"+project_name+"/data/manifest/",exist_ok=True)
    os.makedirs("../../csl_results/"+project_name+"/log/",exist_ok=True)
    
    if os.path.exists("../../csl_results/"+project_name+"/log/output_listFastq.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_listFastq.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_listFastq.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_listFastq.txt")
    if param == "n":
        read_path_destination = "../../csl_results/"+project_name+"/data/fastq/"
        print("Read files will be copied to ../../csl_results/"+project_name+"/data/fastq/")
        print("==================================")
        print("Here's the list of available genomes:")
        genome_list = pd.Series(['human','mouse'])
        print("Index")
        print(genome_list)
        print("==================================")
        print("Specify the index to the genome:(e.g 0)")
        print("\033[91m"+"If you want to use our example dataset, type the number"+"\033[94m"+" 1"+"\x1b[0m")
        genome_index = int(input())
        genome = genome_list[genome_index]
        print("==================================")
        print("Copy the path to your original read files folder:")
        print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
        print("../scripts_DoNotTouch/test/fastq/")
        read_path_original = input()
        read_path_original = os.path.expanduser(read_path_original)
        print("==================================")
        print("Copy the path to design matrix folder (If it's in your home folder, type tilde sign ~):")
        print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
        print("../scripts_DoNotTouch/test/manifest/")
        inpath_design = input()
        inpath_design = os.path.expanduser(inpath_design)
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

        conf = open("../scripts_DoNotTouch/config.py", "w")
        conf.write("project_name="+"'"+project_name+"'"+"\n")
        conf.write("parameters_exist="+"'"+param+"'"+"\n")
        conf.write("read_path_original="+"'"+read_path_original+"'"+"\n")
        conf.write("read_path_destination="+"'"+read_path_destination+"'"+"\n")
        conf.write("genome="+"'"+genome+"'"+"\n")
        conf.write("pairing="+"'"+pairing+"'"+"\n")
        conf.write("inpath_design="+"'"+inpath_design+"'"+"\n")
        conf.write("scriptpath_listdir="+"'"+scriptpath_listdir+"'"+"\n")
        conf.write("scriptpath_copy="+"'"+scriptpath_copy+"'")
        conf.close()

    else:
    
        read_path_original=config.read_path_original
        read_path_destination=config.read_path_destination
        genome=config.genome
        pairing=config.pairing
        inpath_design=config.inpath_design
        scriptpath_listdir=config.scriptpath_listdir
        scriptpath_copy=config.scriptpath_copy
    
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
    print(job[1])
    jobid.append(job[1].split(' ')[2])
    
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
    dirlist.index = range(len(dirlist))
    
    for i in range(len(dirlist)):
        for j in range(len(filename)):
            if dirlist[i] in filename[j]:
                #newname=des.iloc[j,0]+re.sub(r'^.*?_R', '_R', dirlist[i])
                #os.rename(directory+dirlist[i],directory+newname)
                if dirlist.str.contains("_R1_")[i]:
                    newname=des.iloc[j,0]+re.sub(r'^.*?_R1', '_R1', dirlist[i])
                    os.rename(directory+dirlist[i],directory+newname)
                elif dirlist.str.contains("_R2_")[i]:
                    newname=des.iloc[j,0]+re.sub(r'^.*?_R2', '_R2', dirlist[i])
                    os.rename(directory+dirlist[i],directory+newname)
     
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
    read_path_destination = os.path.expanduser(read_path_destination)
    print("========================================")
    
    return read_path_destination+"/"

def fastqc_RunQC(readlist,outdir_fastqc,read_path_destination,scriptpath_fastqc):
            
    global project_name
    
    jobid = []
    for file in readlist:
        command = "source "+scriptpath_fastqc+" "+read_path_destination+file+" "+outdir_fastqc+"/."+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job[1])
        
        jobid.append(job[1].split(' ')[2])
    
    return jobid

def fastqc_ListDir(outdir_fastqc):
    
    dirlist = pd.Series(os.listdir(outdir_fastqc))
    dirlist = dirlist[dirlist.str.endswith('.html')]
    dirlist.index = range(len(dirlist))
    
    DataFrame(dirlist).to_html(outdir_fastqc+'/QC_list.html')
    imgkit.from_file(outdir_fastqc+'/QC_list.html', outdir_fastqc+'/QC_list.png')
    
    #dir_html = IFrame(outdir_fastqc+'/QC_list.html', width=1000, height=800)
    #dir_html = HTML(open(outdir_fastqc+'/QC_list.html').read())
    dir_html = Image(outdir_fastqc+'/QC_list.png')
    os.remove(outdir_fastqc+'/QC_list.png')
    
    return dir_html

def fastqc_Visualization(outdir_fastqc):

    dirlist = pd.Series(os.listdir(outdir_fastqc))
    dirlist = dirlist[dirlist.str.endswith('.html')]
    dirlist.index = range(len(dirlist))
    #print(dirlist)
    
    print("Specify index to visualize HTML file:(e.g 0)")
    index_files = int(input())    
    
    file = dirlist[index_files]
    options = {
    'format': 'png',
    'width': 1500,  # adjust as needed
    'height': 8000  # large enough to capture full content
    }
    imgkit.from_file(outdir_fastqc+file, outdir_fastqc+file+'.png',options=options)
    
    #qc = IFrame(outdir_fastqc+file, width=1000, height=800)
    #qc = HTML(open(outdir_fastqc+file).read())
    qc = Image(outdir_fastqc+file+'.png')
    os.remove(outdir_fastqc+file+'.png')
    
    return qc

def cutadapt_Prep(directory,pairing):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_cutadapt.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_cutadapt.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_cutadapt.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_cutadapt.txt")
    
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
    read_path_destination = os.path.expanduser(read_path_destination)
    print("========================================")
    
    return read_path_destination+"/"

def cutadapt_RunTrimming(adapter,adapter2,minlen,read1_list,read2_list,trimmed1_list,trimmed2_list,scriptpath_cutadapt):
            
    global project_name
    
    jobid = []
    for i in range(len(read1_list)):
        command = "source "+scriptpath_cutadapt+" "+minlen+" "+adapter+" "+adapter2+" "+trimmed1_list[i]+" "+trimmed2_list[i]+" "+read1_list[i]+" "+read2_list[i]+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job[1])
        
        jobid.append(job[1].split(' ')[2])
    
    return jobid

def star_Prep(genome,pairing,read_dir):
        
    global project_name

    if os.path.exists("../../csl_results/"+project_name+"/log/output_star.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_star.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_star.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_star.txt")
    
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        read_dir = "../../csl_results/"+project_name+"/data/cutadapt/"
    print("========================================")
    
    prefix = pd.Series(os.listdir(read_dir))
    prefix = prefix[prefix.str.endswith('fastq.gz')]
    prefix = prefix.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    out_dir = "../../csl_results/"+project_name+"/data/star/"
   
    for i in range(len(prefix)):
        os.makedirs(out_dir+prefix[i],exist_ok=True)

    print("STAR alignment results will be stored in ../../csl_results/"+project_name+"/data/star/")
    
    if genome == 'mouse':
        genome_index_path = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/GRCm39_M29_gencode_starindex"
    elif genome == 'human':
        genome_index_path = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/hg38_p13_gencode_rel42_all_starindex"
    
    read1_list = read_dir+'/'+prefix+'_R1_001.fastq.gz'
    read2_list = read_dir+'/'+prefix+'_R2_001.fastq.gz'
    out_prefix_list = out_dir+prefix+'/'+prefix
    
    if pairing == 'y':
        scriptpath_star = '../scripts_DoNotTouch/STAR/qsub_star_PE.sh'
    else:
        scriptpath_star = '../scripts_DoNotTouch/STAR/qsub_star_SE.sh'

    return genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,scriptpath_star

def star_PrepDirect():
    
    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Are the reads paired-end:(e.g y/n)")
    pairing = input()
    print("========================================")
    print("Specify the path to fastq folder used for alignment:")
    read_path_destination = input()
    read_path_destination = os.path.expanduser(read_path_destination)
    print("========================================")
    
    return genome,pairing,read_path_destination+"/"

def star_RunAlignment(genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,scriptpath_star):
        
    global project_name
    
    jobid = []
    for i in range(len(out_prefix_list)):
        command = "source "+scriptpath_star+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[1])
        jobid.append(job[1].split(' ')[2])
    
    return jobid

def star_ListDir(directory):
    
    global project_name
    starlogdir = "../../csl_results/"+project_name+"/data/star_summary/"
    os.makedirs(starlogdir,exist_ok=True)
    
    dirlist = pd.Series(os.listdir(directory))
    
    i = 0

    for file in dirlist:
        i += 1
        
        logfile = directory+file+"/"+file+"Log.final.out"

        log_df = pd.read_table(logfile,comment='#',header=None,sep="\t",index_col=[0])
        log_raw = log_df.rename({log_df.columns[0]:file},axis='columns')
        if i == 1 :
            log_matrix = log_raw
        else :
            log_matrix = pd.concat([log_matrix,log_raw],axis=1)

    log_matrix.to_csv(starlogdir+'summary_matrix.txt',sep='\t')        
    
    return log_matrix

def kallisto_Prep(genome,pairing,read_dir):
        
    global project_name

    if os.path.exists("../../csl_results/"+project_name+"/log/output_kallisto.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_kallisto.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_kallisto.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_kallisto.txt")
    
    print("========================================")
    print("If you have trimmed the adapters with cutadapt prior, do you want to use these trimmed reads instead?:(y/n)")
    usetrim = input()
    if usetrim == 'y':
        read_dir = "../../csl_results/"+project_name+"/data/cutadapt/"
    print("========================================")
    
    prefix = pd.Series(os.listdir(read_dir))
    prefix = prefix[prefix.str.endswith('fastq.gz')]
    prefix = prefix.str.replace('_R1_001.fastq.gz','',regex=False).str.replace('_R2_001.fastq.gz','',regex=False).unique()
    
    out_dir_kal = "../../csl_results/"+project_name+"/data/kallisto/"
   
    for i in range(len(prefix)):
        os.makedirs(out_dir_kal+prefix[i],exist_ok=True)

    print("Kallisto pseudo-alignment results will be stored in ../../csl_results/"+project_name+"/data/kallisto/")
    
    if genome == 'mouse':
        genome_index_path = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.transcripts.idx"
    elif genome == 'human':
        genome_index_path = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v45.transcripts.idx"
    
    read1_list = read_dir+'/'+prefix+'_R1_001.fastq.gz'
    read2_list = read_dir+'/'+prefix+'_R2_001.fastq.gz'
    out_prefix_list = out_dir_kal+prefix+'/'
    
    if pairing == 'y':
        scriptpath_kallisto = '../scripts_DoNotTouch/Kallisto/qsub_kallisto_PE.sh'
    else:
        scriptpath_kallisto = '../scripts_DoNotTouch/Kallisto/qsub_kallisto_SE.sh'

    return genome_index_path,read1_list,read2_list,out_prefix_list,out_dir_kal,scriptpath_kallisto

def kallisto_PrepDirect():
    
    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Are the reads paired-end:(e.g y/n)")
    pairing = input()
    print("========================================")
    print("Specify the path to fastq folder used for alignment:")
    read_path_destination = input()
    read_path_destination = os.path.expanduser(read_path_destination)
    print("========================================")
    
    return genome,pairing,read_path_destination+"/"

def kallisto_RunAlignment(genome_index_path,read1_list,read2_list,out_prefix_list,out_dir_kal,scriptpath_kallisto):
        
    global project_name
    
    jobid = []
    for i in range(len(out_prefix_list)):
        command = "source "+scriptpath_kallisto+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[1])
        jobid.append(job[1].split(' ')[2])
    
    return jobid

def featurecounts_ListDir(prefix,count_prefix_list):
    
    global project_name
    outpath_counts = "../../csl_results/"+project_name+"/data/counts/"
    os.makedirs(outpath_counts,exist_ok=True)
    print("Featurecounts summary matrix is stored in ../../csl_results/"+project_name+"/data/counts/")
    
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

def featurecounts_Prep(genome,out_dir,pairing):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_featurecounts.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_featurecounts.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_featurecounts.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_featurecounts.txt")
    
    prefix = pd.Series(os.listdir(out_dir))
    
    count_dir = "../../csl_results/"+project_name+"/data/featurecounts/"
   
    for i in range(len(prefix)):
        os.makedirs(count_dir+prefix[i],exist_ok=True)
    
    if genome == 'mouse':
        GTF = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation.gtf"
        strandBED = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation_forStrandDetect_geneID.bed"
    elif genome == 'human':
        GTF = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf"
        strandBED = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_forStrandDetect_geneID.bed"
    
    if pairing == "y":
        scriptpath_featurecounts = '../scripts_DoNotTouch/featureCounts/qsub_featurecounts_PE.sh'
    else:
        scriptpath_featurecounts = '../scripts_DoNotTouch/featureCounts/qsub_featurecounts_SE.sh'

    print("Specify the genomic feature to quantify (e.g gene_name, gene_id, etc):")
    feature = input() 
    
    count_prefix_list = count_dir+prefix+'/'+prefix
    bam_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByCoord.out.bam'

    return scriptpath_featurecounts,GTF,bam_list,count_prefix_list,prefix,feature,strandBED

def featurecounts_PrepDirect():
    
    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Are the reads paired-end:(e.g y/n)")
    pairing = input()
    print("========================================")
    print("Specify the path to alignment folder used for quantification:")
    out_dir = input()
    out_dir = os.path.expanduser(out_dir)
    print("========================================")
    
    return genome,pairing,out_dir+"/"

def featurecounts_RunQuantification(scriptpath_featurecounts,GTF,bam_list,count_prefix_list,feature,strandBED):
     
    global project_name
    
    jobid = []
    for i in range(len(bam_list)):
        command = "source "+scriptpath_featurecounts+" "+bam_list[i]+" "+GTF+" "+feature+" "+count_prefix_list[i]+" "+strandBED+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[1])
        jobid.append(job[1].split(' ')[2])
    
    return jobid

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

def rsem_Prep(genome,out_dir,pairing):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_rsem.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_rsem.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_rsem.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_rsem.txt")
    
    prefix = pd.Series(os.listdir(out_dir))
    
    count_dir = "../../csl_results/"+project_name+"/data/rsem/"
   
    for i in range(len(prefix)):
        os.makedirs(count_dir+prefix[i],exist_ok=True)
    
    if genome == 'mouse':
        rsem_index = "/grid/bsr/data/data/utama/genome/mouse_rsem_index_star_gencode_GRCm39_M29_v2.7.10a/mouse_gencode"
        strandBED = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation_forStrandDetect_geneID.bed"
    elif genome == 'human':
        rsem_index = "/grid/bsr/data/data/utama/genome/human_rsem_index_star_gencode_hg38_p13_rel42_v2.7.2b/human_gencode"
        strandBED = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_forStrandDetect_geneID.bed"
    
    if pairing == "y":
        scriptpath_rsem = '../scripts_DoNotTouch/RSEM/qsub_RSEM_PE.sh'
    else:
        scriptpath_rsem = '../scripts_DoNotTouch/RSEM/qsub_RSEM_SE.sh'

    feature = "gene_id"
    
    count_prefix_list = count_dir+prefix+'/'+prefix
    bam_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByCoord.out.bam'
    bamTranscript_list = out_dir+prefix+'/'+prefix+'Aligned.toTranscriptome.out.bam'

    return scriptpath_rsem,rsem_index,bam_list,count_prefix_list,prefix,feature,strandBED,bamTranscript_list

def rsem_PrepDirect():
    
    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Are the reads paired-end:(e.g y/n)")
    pairing = input()
    print("========================================")
    print("Specify the path to alignment folder used for quantification:")
    out_dir = input()
    out_dir = os.path.expanduser(out_dir)
    print("========================================")
    
    return genome,pairing,out_dir+"/"

def rsem_RunQuantification(scriptpath_rsem,rsem_index,bam_list,count_prefix_list,feature,strandBED,bamTranscript_list):
     
    global project_name
    
    jobid = []
    for i in range(len(bam_list)):
        command = "source "+scriptpath_rsem+" "+bam_list[i]+" "+rsem_index+" "+feature+" "+count_prefix_list[i]+" "+strandBED+" "+bamTranscript_list[i]+" "+project_name
        job = os.popen(command).read().splitlines()
        print(job[1])
        jobid.append(job[1].split(' ')[2])
    
    return jobid
    
def rsem_CreateCountMatrix():

    global project_name
    inpath_counts = "../../csl_results/"+project_name+"/data/rsem/"
    outpath_counts = "../../csl_results/"+project_name+"/data/counts/"
    print("Count matrix is stored in ../../csl_results/"+project_name+"/data/counts/")

    filelist=sorted([f for f in os.listdir(inpath_counts) if not f.startswith('.')])

    i = 0

    for file in filelist:
        i += 1
        gene_df = pd.read_table(os.path.join(inpath_counts,file,file+'.genes.results'),comment='#',header=[0],index_col=[0])
        isoform_df = pd.read_table(os.path.join(inpath_counts,file,file+'.isoforms.results'),comment='#',header=[0],index_col=[0])
        
        gene_tpm_raw = gene_df.drop(['transcript_id(s)','length','effective_length','expected_count','FPKM'],axis=1)
        gene_fpkm_raw = gene_df.drop(['transcript_id(s)','length','effective_length','expected_count','TPM'],axis=1)
        
        isoform_tpm_raw = isoform_df.drop(['length','effective_length','expected_count','IsoPct','FPKM'],axis=1)
        isoform_fpkm_raw = isoform_df.drop(['length','effective_length','expected_count','IsoPct','TPM'],axis=1)
        
        gene_tpm_raw = gene_tpm_raw.rename({gene_tpm_raw.columns[0]:file},axis='columns')
        gene_fpkm_raw = gene_fpkm_raw.rename({gene_fpkm_raw.columns[0]:file},axis='columns')
        isoform_tpm_raw = isoform_tpm_raw.rename({isoform_tpm_raw.columns[1]:file},axis='columns')
        isoform_fpkm_raw = isoform_fpkm_raw.rename({isoform_fpkm_raw.columns[1]:file},axis='columns')
        
        if i == 1 :
            gene_tpm_matrix = gene_tpm_raw
            gene_fpkm_matrix = gene_fpkm_raw
            isoform_tpm_matrix = isoform_tpm_raw
            isoform_fpkm_matrix = isoform_fpkm_raw
        else :
            gene_tpm_matrix = pd.concat([gene_tpm_matrix,gene_tpm_raw],axis=1)
            gene_fpkm_matrix = pd.concat([gene_fpkm_matrix,gene_fpkm_raw],axis=1)
            
            isoform_tpm_raw = isoform_tpm_raw.drop(['gene_id'],axis=1)
            isoform_fpkm_raw = isoform_fpkm_raw.drop(['gene_id'],axis=1)
            
            isoform_tpm_matrix = pd.concat([isoform_tpm_matrix,isoform_tpm_raw],axis=1)
            isoform_fpkm_matrix = pd.concat([isoform_fpkm_matrix,isoform_fpkm_raw],axis=1)

    gene_tpm_matrix.to_csv(outpath_counts+'gene_tpm_matrix.txt',sep='\t')
    gene_fpkm_matrix.to_csv(outpath_counts+'gene_fpkm_matrix.txt',sep='\t')
    
    isoform_tpm_matrix.to_csv(outpath_counts+'isoform_tpm_matrix.txt',sep='\t')
    isoform_fpkm_matrix.to_csv(outpath_counts+'isoform_fpkm_matrix.txt',sep='\t')

    return outpath_counts,gene_tpm_matrix,gene_fpkm_matrix,isoform_tpm_matrix,isoform_fpkm_matrix

def deseq2_Prep(inpath_design):
    
    global project_name
    outpath = "../../csl_results/"+project_name+"/data/deseq2/"
    os.makedirs(outpath,exist_ok=True)
    print("DESeq2 results are stored in ../../csl_results/"+project_name+"/data/deseq2/")

    design = pd.read_table(inpath_design+'/design_matrix.txt',index_col=0)
    design = design.iloc[:,:len(design.columns)-1]

    print("========================================")
    print("If you have a redundant variable/column in your design matrix,")
    print("please type the name of that variable/column(e.g age)")
    print("Otherwise, leave it blank and press enter/return")
    redundant=input()
    if len(redundant) == 0:
        redundant = "NoRedundant"
    
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
    
    return scriptpath_deseq2,Rpath_deseq2,outpath,refcond,compared,design_var,redundant

def deseq2_PrepDirect():
    
    print("========================================")
    print("Specify the path to folder containing count_matrix.txt used for DE:")
    outpath_counts = input()
    outpath_counts = os.path.expanduser(outpath_counts)
    print("========================================")
    print("Specify the path to folder containing design_matrix.txt used for DE:")
    inpath_design = input()
    inpath_design = os.path.expanduser(inpath_design)
    print("========================================")
    
    return outpath_counts+"/",inpath_design+"/"

def deseq2_RunDE(scriptpath_deseq2,Rpath_deseq2,inpath_counts,inpath_design,outpath,refcond,compared,redundant):
    
    global project_name
    
    jobid = []
    command = "source "+scriptpath_deseq2+" "+Rpath_deseq2+" "+inpath_counts+"/count_matrix.txt"+" "+inpath_design+"/design_matrix.txt"+" "+outpath+" "+refcond+" "+compared+" "+redundant+" "+project_name
    job = os.popen(command).read().splitlines()
    print(job[1])
    jobid.append(job[1].split(' ')[2])

    return jobid

def gseapy_Prep():
    
    global project_name
    
    outpath_pathway = "../../csl_results/"+project_name+"/data/gseapy/"
    os.makedirs(outpath_pathway,exist_ok=True)
    print("GSEApy results are stored in ../../csl_results/"+project_name+"/data/gseapy/")

    print("========================================")
    print("Specify gene set database:")
    print("MSigDB_Hallmark_2020, KEGG_2021_Human, GO_Biological_Process_2025, Reactome_Pathways_2024")
    print("ARCHS4_TFs_Coexp, ENCODE_TF_ChIP-seq_2015, ENCODE_Histone_Modifications_2015")
    print("FANTOM6_lncRNA_KD_DEGs, miRTarBase_2017, TRANSFAC_and_JASPAR_PWMs")
    print("GTEx_Tissues_V8_2023, CellMarker_2024, Cancer_Cell_Line_Encyclopedia")
    print("ClinVar_2019, GTEx_Aging_Signatures_2021, Proteomics_Drug_Atlas_2023")
    geneset = input()
    
    return geneset,outpath_pathway

def gseapy_PrepDirect():

    global project_name

    outpath_pathway = "../../csl_results/"+project_name+"/data/gseapy/"
    os.makedirs(outpath_pathway,exist_ok=True)
    print("GSEApy results are stored in ../../csl_results/"+project_name+"/data/gseapy/")

    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    
    print("========================================")
    print("Specify the genomic feature to quantify (e.g gene_name, gene_id, etc):")
    feature = input() 
    
    print("========================================")
    print("Specify the path to folder containing design_matrix.txt used for DE:")
    inpath_design = input()
    inpath_design = os.path.expanduser(inpath_design)

    print("========================================")
    print("Specify the path to folder containing normalized counts from DE:")
    outpath = input()

    print("========================================")
    print("Which phenotype/condition/replicate/batch should be the reference/baseline?(e.g control)")
    refcond = input()
    
    print("========================================")
    print("Which phenotype/condition/replicate/batch to compare?(e.g treated)")
    compared = input()

    print("========================================")
    print("Specify gene set database:")
    print("MSigDB_Hallmark_2020, KEGG_2021_Human, GO_Biological_Process_2025, Reactome_Pathways_2024")
    print("ARCHS4_TFs_Coexp, ENCODE_TF_ChIP-seq_2015, ENCODE_Histone_Modifications_2015")
    print("FANTOM6_lncRNA_KD_DEGs, miRTarBase_2017, TRANSFAC_and_JASPAR_PWMs")
    print("GTEx_Tissues_V8_2023, CellMarker_2024, Cancer_Cell_Line_Encyclopedia")
    print("ClinVar_2019, GTEx_Aging_Signatures_2021, Proteomics_Drug_Atlas_2023")
    geneset = input()
    
    return geneset,genome,feature,inpath_design+"/",outpath+"/",outpath_pathway,refcond,compared

def gseapy_RunPathway(geneset,genome,feature,inpath_design,outpath,outpath_pathway,refcond,compared):
    
    global project_name
    
    design = pd.read_table(inpath_design+'/design_matrix.txt',index_col=0)
    vardesign = design.T[(design.T.iloc[:,0]==refcond) | (design.T.iloc[:,0]==compared)].T.columns[0]
    design = design.loc[design[vardesign].isin([refcond,compared]),:]
    class_vector = list(design[vardesign])
    reordering = ['GENE','NAME']+list(design.index)

    gene_exp = pd.read_table(outpath+'/normalized_counts_'+compared+'_vs_'+refcond+'(ref).txt',index_col=0)
    gene_exp = gene_exp.drop(['DESCRIPTION'],axis=1)
    gene_exp.index = gene_exp.index.str.split('.').str[0]
    gene_exp['GENE'] = gene_exp.index
    gene_exp['NAME'] = 'na'

    #######################

    bm = Biomart()
    # note the dataset and attribute names are different
    m2h = bm.query(dataset='mmusculus_gene_ensembl',
                   attributes=['ensembl_gene_id','external_gene_name',
                               'hsapiens_homolog_ensembl_gene',
                               'hsapiens_homolog_associated_gene_name'])
    h2m = bm.query(dataset='hsapiens_gene_ensembl',
                   attributes=['ensembl_gene_id','external_gene_name',
                               'mmusculus_homolog_ensembl_gene',
                               'mmusculus_homolog_associated_gene_name'])
    
    if (genome == 'mouse') & (feature == 'gene_name'):        
        gene_exp_conv = gene_exp.merge(m2h,how='inner',left_index=True,right_on='external_gene_name')
        gene_exp_conv['GENE'] = gene_exp_conv['hsapiens_homolog_associated_gene_name']
    elif (genome == 'mouse') & (feature == 'gene_id'):
        gene_exp_conv = gene_exp.merge(m2h,how='inner',left_index=True,right_on='ensembl_gene_id')
        gene_exp_conv['GENE'] = gene_exp_conv['hsapiens_homolog_associated_gene_name']       
    elif (genome == 'human') & (feature == 'gene_id'):
        gene_exp_conv = gene_exp.merge(h2m,how='inner',left_index=True,right_on='ensembl_gene_id')
        gene_exp_conv['GENE'] = gene_exp_conv['external_gene_name']
    elif (genome == 'human') & (feature == 'gene_name'):
        gene_exp_conv = gene_exp
    
    gene_exp_conv = gene_exp_conv.dropna(subset=['GENE'])
    gene_exp_conv = gene_exp_conv[reordering]

    #############################
    
    gs = GSEA(data=gene_exp_conv,
         gene_sets=geneset,
         classes = class_vector, # cls=class_vector
         # set permutation_type to phenotype if samples >=15
         permutation_type='gene_set',
         permutation_num=1000, # reduce number to speed up test
         outdir=outpath_pathway,
         method='signal_to_noise',
         threads=4, seed= 8)
    
    gs.pheno_pos = compared
    gs.pheno_neg = refcond
    gs_res = gs.run()

    pathways = pd.read_csv(outpath_pathway+'/gseapy.gene_set.gsea.report.csv',index_col=None)
    pathways = pathways.sort_values(by=['FDR q-val'])
    pathways.to_csv(outpath_pathway+'/report.gseapy.'+geneset+'.csv',index=None)
    terms = pathways.Term

    return gs,gs_res,pathways,terms,project_name

def gseapy_DotPlot(outpath_pathway,pathways,geneset):
    
    # to save your figure, make sure that ``ofname`` is not None
    dot = dotplot(pathways,
             column="FDR q-val",
             title=geneset,
             cmap=plt.cm.YlOrRd,
             size=8,
             ofname=outpath_pathway+"DotPlot_Top10."+geneset+'.png',
             figsize=(4,5), cutoff=1)
    dot = dotplot(pathways,
             column="FDR q-val",
             title=geneset,
             cmap=plt.cm.YlOrRd,
             size=8,
             figsize=(4,5), cutoff=1)
    
    return dot

def gseapy_EnrichPlot(pathways,gs):

    terms = pathways.Term

    print("Specify index to visualize Enrichment plot of a selected pathway:(e.g 0)")
    index_files = int(input())

    enrich = gs.plot(terms[index_files]) # If choosing only 1 pathway

    return enrich

def gseapy_heatmap(pathways,gs):

    terms = pathways.Term
    
    print("Specify index to visualize Heatmap plot of a selected pathway:(e.g 0)")
    index_files = int(input())

    genes = pathways.Lead_genes[index_files].split(";")
    # Make sure that ``ofname`` is not None, if you want to save your figure to disk
    heatgsea = heatmap(df = gs.heatmat.loc[genes], z_score=0, title=terms[index_files], figsize=(14,4))

    return heatgsea

def visualization_PrepDirect():
    
    print("========================================")
    print("Specify the path to folder containing design_matrix.txt used for DE:")
    inpath_design = input()
    inpath_design = os.path.expanduser(inpath_design)
    print("========================================")
    print("Specify the path to folder containing DESeq2 results:")
    outpath = input()
    outpath = os.path.expanduser(outpath)
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
    count_norm = pd.read_table(outpath+'/normalized_counts_'+compared+'_vs_'+refcond+'(ref).txt',index_col=0)
    count_norm = count_norm.drop(['DESCRIPTION'],axis=1)
    count_norm_sig = count_norm[count_norm.index.isin(genes)]
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300 
    for i in range(len(design.columns)):
        lut = dict(zip(design.iloc[:,i].unique(), sns.color_palette("Pastel2")))
        col_colors = design.iloc[:,i].map(lut)
        heat = sns.clustermap(count_norm_sig,z_score=0,cmap='vlag',col_colors=col_colors)
        heat.savefig(outpath+"/heatmap_"+design.columns[i]+"_"+compared+'_vs_'+refcond+'(ref).png')

    return heat

def visualization_pca(design_var):

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
