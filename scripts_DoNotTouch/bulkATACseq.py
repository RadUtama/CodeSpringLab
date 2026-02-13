
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
import imgkit

project_name=config.project_name
param=config.parameters_exist

def Tree():
    
    cmd = "tree ../../csl_results/"
    print(os.popen(cmd).read())

def DeleteJobs(jobid):
    
    for i in range(len(jobid)):
        command = "scancel "+jobid[i]
        #command = "qdel "+jobid[i]
        jobdel = os.popen(command).read().splitlines()
        print('Deleting job ID: '+jobid[i])
        while True:
            command = "squeue -u $USER"
            #command = "qstat"
            del_status = os.popen(command).read().splitlines()
            if len(del_status) == 1:
            #if del_status == []:
                break
            else:
                del_stat = pd.DataFrame(del_status)[0][1:].str.split(expand=True)[0].str.contains(jobid[i]).any()
                if del_stat == False:
                    break
            time.sleep(5)
            
def DeleteOneJob(jobid):

    command = "scancel "+jobid[0]
    #command = "qdel "+jobid[0]
    jobdel=os.popen(command).read().splitlines()
    print('Deleting job ID: '+jobid[0])
    while True:
        command = "squeue -u $USER"
        #command = "qstat"
        del_status=os.popen(command).read().splitlines()
        if len(del_status) == 1:
        #if del_status == []:
            break
        else:
            del_stat = pd.DataFrame(del_status)[0][1:].str.split(expand=True)[0].str.contains(jobid[0]).any()
            if del_stat == False:
                break
        time.sleep(5)
            
def Qstat(jobid):
        
    print("==================================")
    starttime = time.time()
    while True:
        command = "squeue -u $USER"
        #command = "qstat"
        running_status = os.popen(command).read().splitlines()
    
        clear_output(wait=True)

        if len(running_status) == 1:
        #if running_status == []:
            print("All jobs done !" )
            print(".........................")
            print("Running time: "+str( round( ((time.time() - starttime) / 60) , 2 ) )+" minutes")
            runstat = "done"
        elif pd.DataFrame(running_status)[0][1:].str.split(expand=True)[0].str.contains('|'.join(jobid)).any():
        #elif pd.DataFrame(running_status)[0][2:].str.split(expand=True)[0].str.contains('|'.join(jobid)).any():

            all_stat = pd.DataFrame(running_status)[0][1:].str.split(expand=True)[4]            
            #all_stat = pd.DataFrame(running_status)[0][2:].str.split(expand=True)[4]
            id_sub = pd.DataFrame(running_status)[0][1:].str.split(expand=True)[0].str.contains('|'.join(jobid))
            #id_sub = pd.DataFrame(running_status)[0][2:].str.split(expand=True)[0].str.contains('|'.join(jobid))
            id_stat = all_stat[id_sub]
            status_stat = id_stat.value_counts().rename(index={'R':'Number of jobs still running = ',
                                                'PD':'Number of jobs waiting in line = ',
                                                'CG' : 'Number of jobs still finishing = ',
                                                'CD' : 'Number of jobs completed = ',
                                                'F' : 'Number of jobs failed = ',
                                                'CA' : 'Number of jobs cancelled = ',
                                                'TO' : 'Number of jobs timeout = ',
                                                'CF' : 'Number of jobs configuring = '
                                                })
            #status_stat = id_stat.value_counts().rename(index={'r':'Number of jobs still running = ',
            #                                        'qw':'Number of jobs waiting in line = ',
            #                                        't':'Number of jobs about to run = ',
            #                                        'Rt':'Number of jobs about to run = ',
            #                                        'Rr':'Number of jobs about to run = ',
            #                                        'Eqw':'Number of jobs cannot run (server error) = '
            #                                       })
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
        print("\033[91m"+"If you want to use our example dataset, type the number"+"\033[94m"+" 0"+"\x1b[0m")
        genome_index = int(input())
        genome = genome_list[genome_index]
        print("==================================")
        print("Copy the path to your original read files folder:")
        print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
        print("../scripts_DoNotTouch/test/fastq_atac/")
        read_path_original = input()
        read_path_original = os.path.expanduser(read_path_original)
        print("==================================")
        print("Copy the path to design matrix folder (If it's in your home folder, type tilde sign ~):")
        print("\033[91m"+"If you want to use our example dataset, copy and paste this path below,"+"\x1b[0m")
        print("../scripts_DoNotTouch/test/manifest_atac/")
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
    
    #command = "sbatch "+scriptpath_listdir+" "+read_path_original+" "+project_name
    ##command = "source "+scriptpath_listdir+" "+read_path_original+" "+project_name
    ##joblist=os.popen(command).read().splitlines()
    ##print(joblist)
    
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
    read_path_destination = os.path.expanduser(read_path_destination)
    print("========================================")
    
    return read_path_destination+"/"

def filetransfer_Copy(read_path_original,scriptpath_copy):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/data/fastq") :
        shutil.rmtree("../../csl_results/"+project_name+"/data/fastq")
    
    jobid = []

    stderr = "-e ../../csl_results/"+project_name+"/log/error_copyFastq.txt"
    stdout = "-o ../../csl_results/"+project_name+"/log/output_copyFastq.txt"
    command = "sbatch "+stderr+" "+stdout+" "+scriptpath_copy+" "+read_path_original+" "+"../../csl_results/"+project_name+"/data/fastq"+" "+project_name
    #command = "source "+scriptpath_copy+" "+read_path_original+" "+"../../csl_results/"+project_name+"/data/fastq"+" "+project_name
    job = os.popen(command).read().splitlines()
    print(job[0])
    #print(job[1])
    jobid.append(job[0].split(' ')[3])
    #jobid.append(job[1].split(' ')[2])
    
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
            if dirlist.iloc[i] in filename.iloc[j]:
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
        stderr = "-e ../../csl_results/"+project_name+"/log/error_fastQC.txt"
        stdout = "-o ../../csl_results/"+project_name+"/log/output_fastQC.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_fastqc+" "+read_path_destination+file+" "+outdir_fastqc+"/."+" "+project_name
        #command = "source "+scriptpath_fastqc+" "+read_path_destination+file+" "+outdir_fastqc+"/."+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])

        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
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
        stderr = "-e ../../csl_results/"+project_name+"/log/error_cutadapt.txt"
        stdout = "-o ../../csl_results/"+project_name+"/log/output_cutadapt.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_cutadapt+" "+minlen+" "+adapter+" "+adapter2+" "+trimmed1_list[i]+" "+trimmed2_list[i]+" "+read1_list[i]+" "+read2_list[i]+" "+project_name
        #command = "source "+scriptpath_cutadapt+" "+minlen+" "+adapter+" "+adapter2+" "+trimmed1_list[i]+" "+trimmed2_list[i]+" "+read1_list[i]+" "+read2_list[i]+" "+project_name 
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])

        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid

def bowtie2_Prep(genome,pairing,read_dir):
        
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_bowtie2.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_bowtie2.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_bowtie2.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_bowtie2.txt")
    
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
        effgenomesize = "2654621783"
        chromsize = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/chrom.sizes"
        genome_index_path = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/GRCm39_M29_gencode_bowtie2index/GRCm39_M29_gencode"
    elif genome == 'human':
        effgenomesize = "2913022398"
        chromsize = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/chrom.sizes"
        genome_index_path = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/hg38_p13_gencode_bowtie2index/hg38_p13_gencode"
    
    read1_list = read_dir+'/'+prefix+'_R1_001.fastq.gz'
    read2_list = read_dir+'/'+prefix+'_R2_001.fastq.gz'
    out_prefix_list = out_dir+prefix+'/'+prefix
    
    if pairing == 'y':
        scriptpath_bowtie2 = '../scripts_DoNotTouch/bowtie2/qsub_bowtie2_PE.sh'
    else:
        scriptpath_bowtie2 = '../scripts_DoNotTouch/bowtie2/qsub_bowtie2_SE.sh'

    return genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,effgenomesize,chromsize,scriptpath_bowtie2

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
    read_path_destination = os.path.expanduser(read_path_destination)
    print("========================================")
    
    return genome,pairing,read_path_destination+"/"

def bowtie2_RunAlignment(genome_index_path,read1_list,read2_list,out_prefix_list,out_dir,effgenomesize,chromsize,scriptpath_bowtie2):
        
    global project_name
    
    cmd_dir = "cat > ../../csl_results/"+project_name+"/log/error_bowtie2.txt"
    cmd_dir_run = os.popen(cmd_dir)
    
    jobid = []
    for i in range(len(out_prefix_list)):
        stderr = "-e "+out_prefix_list[i]+"Log.final.out"
        stdout = "-o ../../csl_results/"+project_name+"/log/output_bowtie2.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_bowtie2+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+effgenomesize+" "+chromsize+" "+project_name
        #command = "source "+scriptpath_bowtie2+" "+out_prefix_list[i]+" "+genome_index_path+" "+read1_list[i]+" "+read2_list[i]+" "+effgenomesize+" "+chromsize+" "+project_name
        #job = os.popen(command).read().strip().splitlines()
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])

        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
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

        log_df = pd.read_csv(logfile,skiprows=0,nrows=15)
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
    #qc = IFrame(out_dir+"/"+file+"/"+file+"_insert_size_histogram.jpg", width=1000, height=800)
    #qc = HTML(open(out_dir+"/"+file+"/"+file+"_insert_size_histogram.jpg").read())
    qc = Image(out_dir+"/"+file+"/"+file+"_insert_size_histogram.jpg")
    
    return qc

def macs2_Prep(genome,out_dir,pairing):
    
    global project_name
    if os.path.exists("../../csl_results/"+project_name+"/log/output_macs2.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_macs2.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_macs2.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_macs2.txt")
    
    prefix = pd.Series(os.listdir(os.path.expanduser(out_dir)))
    
    macs2_dir = "../../csl_results/"+project_name+"/data/macs2/"
   
    for i in range(len(prefix)):
        os.makedirs(macs2_dir+prefix[i],exist_ok=True)

    print("========================================")
    print("Please set the q-value cutoff (default 0.05):")
    qval = input()
    print("========================================")
    
    if genome == 'mouse':
        homerspecies = "mm39"
        genomesize = "1.87e+9"
        chromsize = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/chrom.sizes"
        anno_onlyChrNoMito = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation_onlyChrNoMito.bed"
    elif genome == 'human':
        homerspecies = "hg38"
        genomesize = "2.7e+9"
        chromsize = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/chrom.sizes"
        anno_onlyChrNoMito = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_onlyChrNoMito.bed"
    
    if pairing == "y":
        scriptpath_macs2 = '../scripts_DoNotTouch/MACS2/qsub_macs2_PE.sh'
    else:
        scriptpath_macs2 = '../scripts_DoNotTouch/MACS2/qsub_macs2_SE.sh'

    print("========================================")
    print("Remove all kinds of read duplicates:(e.g y/n)")
    removeDup = input()
    print("========================================")
    if removeDup == "y":
        bed_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByCoord_removeDup.out.bed'
    else:
        bed_list = out_dir+prefix+'/'+prefix+'Aligned.sortedByCoord.out.bed'
    
    macs2_prefix_list = macs2_dir+prefix+'/'

    return scriptpath_macs2,genomesize,chromsize,bed_list,macs2_prefix_list,prefix,anno_onlyChrNoMito,macs2_dir,qval,homerspecies

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
    out_dir = os.path.expanduser(out_dir)
    print("========================================")
    
    return genome,pairing,out_dir+"/"

def macs2_RunPeakCalling(scriptpath_macs2,genomesize,chromsize,bed_list,macs2_prefix_list,prefix,anno_onlyChrNoMito,qval,homerspecies,out_dir):
     
    global project_name
    
    jobid = []
    for i in range(len(bed_list)):
        stderr = "-e ../../csl_results/"+project_name+"/log/error_macs2.txt"
        stdout = "-o ../../csl_results/"+project_name+"/log/output_macs2.txt"
        command = "sbatch "+stderr+" "+stdout+" "+scriptpath_macs2+" "+prefix[i]+" "+bed_list[i]+" "+genomesize+" "+chromsize+" "+macs2_prefix_list[i]+" "+anno_onlyChrNoMito+" "+qval+" "+project_name+" "+homerspecies+" "+out_dir
        #command = "source "+scriptpath_macs2+" "+prefix[i]+" "+bed_list[i]+" "+genomesize+" "+chromsize+" "+macs2_prefix_list[i]+" "+anno_onlyChrNoMito+" "+qval+" "+project_name+" "+homerspecies+" "+out_dir
        #job = os.popen(command).read().strip().splitlines()
        job = os.popen(command).read().splitlines()
        print(job[0])
        #print(job[1])

        jobid.append(job[0].split(' ')[3])
        #jobid.append(job[1].split(' ')[2])
    
    return jobid

def macs2_PeakList():
    
    global project_name
    
    outpath = "../../csl_results/"+project_name+"/data/macs2/"
    dirlist = DataFrame(pd.Series(os.listdir(outpath)))
    dirlist.index = range(len(dirlist))
    print(dirlist)
    
    print("========================================")
    print("Specify index to visualize insert size frequency:(e.g 0)")
    index_files = int(input())    
    
    file = dirlist.loc[index_files,0]
    
    #peaklist = pd.read_table(outpath+"/"+file+"/"+file+"_peaks.xls",header=0,index_col=0,comment='#',engine='python')
    peaklist = pd.read_table(outpath+"/"+file+"/"+file+"_peaks_annotated.txt",header=0,index_col=0,comment='#',engine='python')
    
    return peaklist

def diffbind_Prep(inpath_design):

    global project_name
    outpath = "../../csl_results/"+project_name+"/data/diffbind/"
    os.makedirs(outpath,exist_ok=True)
    print("DiffBind results are stored in ../../csl_results/"+project_name+"/data/diffbind/")

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
    print("Which phenotype/condition/replicate/batch should be the reference/baseline?(e.g day1)")
    refcond = input()

    print("========================================")
    print("Which phenotype/condition/replicate/batch to compare?(e.g day3)")
    compared = input()

    scriptpath_diffbind = '../scripts_DoNotTouch/DiffBind/qsub_diffbind.sh'
    Rpath_diffbind = '../scripts_DoNotTouch/DiffBind/DiffBind.R'

    return scriptpath_diffbind,Rpath_diffbind,outpath,refcond,compared

def diffbind_PrepDirect():

    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Specify the path to alignment folder used for Differential Peaks:")
    out_dir = input()
    out_dir = os.path.expanduser(out_dir)
    print("========================================")
    print("Specify the path to folder containing .narrowPeak used for Differential Peaks:")
    outpath_peak = input()
    outpath_peak = os.path.expanduser(outpath_peak)
    print("========================================")
    print("Specify the path to folder containing design_matrix.txt used for Differential Peaks:")
    inpath_design = input()
    inpath_design = os.path.expanduser(inpath_design)
    print("========================================")

    return genome,out_dir+"/",outpath_peak+"/",inpath_design+"/"

def diffbind_RunDiffPeak(scriptpath_diffbind,Rpath_diffbind,outpath_peak,inpath_design,outpath,refcond,compared,genome,out_dir):

    global project_name

    bamdir_example = "../../csl_results/example_dataset/data/bowtie2/"
    peakdir_example = "../../csl_results/example_dataset/data/macs2/"
    bamsuffix = "Aligned.sortedByCoord_removeDup.out.bam"
    peaksuffix = "_peaks.narrowPeak"
    
    day1_rep2_bam = bamdir_example+"/CD4_day1_rep2/"
    day3_rep2_bam = bamdir_example+"/CD4_day3_rep2/"
    day1_rep2_peak = peakdir_example+"/CD4_day1_rep2/"
    day3_rep2_peak = peakdir_example+"/CD4_day3_rep2/"
    
    day1_rep3_bam = bamdir_example+"/CD4_day1_rep3/"
    day3_rep3_bam = bamdir_example+"/CD4_day3_rep3/"
    day1_rep3_peak = peakdir_example+"/CD4_day1_rep3/"
    day3_rep3_peak = peakdir_example+"/CD4_day3_rep3/"
    
    os.makedirs(day1_rep3_bam,exist_ok=True)
    os.makedirs(day3_rep3_bam,exist_ok=True)
    os.makedirs(day1_rep3_peak,exist_ok=True)
    os.makedirs(day3_rep3_peak,exist_ok=True)

    shutil.copy(day1_rep2_bam+"/CD4_day1_rep2"+bamsuffix,day1_rep3_bam+"/CD4_day1_rep3"+bamsuffix)
    shutil.copy(day3_rep2_bam+"/CD4_day3_rep2"+bamsuffix,day3_rep3_bam+"/CD4_day3_rep3"+bamsuffix)
    shutil.copy(day1_rep2_peak+"/CD4_day1_rep2"+peaksuffix,day1_rep3_peak+"/CD4_day1_rep3"+peaksuffix)
    shutil.copy(day3_rep2_peak+"/CD4_day3_rep2"+peaksuffix,day3_rep3_peak+"/CD4_day3_rep3"+peaksuffix)
    
    jobid = []
    stderr = "-e ../../csl_results/"+project_name+"/log/error_diffbind.txt"
    stdout = "-o ../../csl_results/"+project_name+"/log/output_diffbind.txt"
    command = "sbatch "+stderr+" "+stdout+" "+scriptpath_diffbind+" "+Rpath_diffbind+" "+outpath+" "+inpath_design+" "+outpath_peak+" "+refcond+" "+compared+" "+genome+" "+out_dir+" "+project_name
    #command = "source "+scriptpath_diffbind+" "+Rpath_diffbind+" "+outpath+" "+inpath_design+" "+outpath_peak+" "+refcond+" "+compared+" "+genome+" "+out_dir+" "+project_name
    job = os.popen(command).read().splitlines()
    print(job[0])
    #print(job[1])

    jobid.append(job[0].split(' ')[3])
    #jobid.append(job[1].split(' ')[2])

    return jobid

def diffbind_PlotCorrHeatByCounts(outpath):

    corrheatbycounts = Image(outpath+"/diffbind_CorrHeatmap_byCounts.png")
    
    return corrheatbycounts

def diffbind_PlotCorrHeatByPeaks(outpath):

    corrheatbypeaks = Image(outpath+"/diffbind_CorrHeatmap_byPeaks.png")
    
    return corrheatbypeaks

def diffbind_PlotSampleHeatByDiffPeaks(outpath):

    sampleheatbydiffpeaks = Image(outpath+"/diffbind_SampleHeatmap_byDiffPeaks.png")
    
    return sampleheatbydiffpeaks

def diffbind_PlotPCAByNormCounts(outpath):

    pcabynormcounts = Image(outpath+"/diffbind_pca_byNormCounts.png")
    
    return pcabynormcounts

def diffbind_PlotPCAByDiffPeaks(outpath):

    pcabydiffpeaks = Image(outpath+"/diffbind_pca_byDiffPeaks.png")
    
    return pcabydiffpeaks

def diffbind_VolcanoByDiffPeaks(outpath):

    volcanobydiffpeaks = Image(outpath+"/diffbind_volcano_byDiffPeaks.png")
    
    return volcanobydiffpeaks

def diffbind_PlotBoxByDiffPeaks(outpath):

    boxbydiffpeaks = Image(outpath+"/diffbind_boxplot_byDiffPeaks.png")
    
    return boxbydiffpeaks

def diffbind_PlotProfileHeatBySampleContrast(outpath):

    profheatbysamplecontrast = Image(outpath+"/diffbind_ProfileHeatmap_bySampleContrast.png")
    
    return profheatbysamplecontrast

def diffbind_PlotProfileHeatByMergedContrast(outpath):

    profheatbymergedcontrast = Image(outpath+"/diffbind_ProfileHeatmap_byMergedContrast.png")
    
    return profheatbymergedcontrast
    
def visualization_PrepDirect():
    
    print("========================================")
    print("Specify genome:(e.g human, mouse, etc)")
    genome = input()
    print("========================================")
    print("Specify the path to MACS2/BED folder from peak calling:")
    out_peak = input()
    out_peak = os.path.expanduser(out_peak)
    print("========================================")
    
    return genome,out_peak+"/"
    
def visualization_Prep(genome,out_peak):
    
    if os.path.exists("../../csl_results/"+project_name+"/log/output_genomeTracks.txt") & os.path.exists("../../csl_results/"+project_name+"/log/error_genomeTracks.txt"):
        os.remove("../../csl_results/"+project_name+"/log/output_genomeTracks.txt")
        os.remove("../../csl_results/"+project_name+"/log/error_genomeTracks.txt")
    
    print("========================================")
    print("Specify which genome region to plot tracks:(e.g chr1:2700000-2800000)")
    region = input()
    print("========================================")
    
    prefix = pd.Series(os.listdir(os.path.expanduser(out_peak)))

    macs2_prefix_list = out_peak+'/'+prefix+'/'
    
    tracks_dir = "../../csl_results/"+project_name+"/data/genomeTracks/"
    os.makedirs(tracks_dir,exist_ok=True)
    
    scriptpath_tracks = '../scripts_DoNotTouch/genomeTracks/qsub_genomeTracks.sh'
    
    if genome == 'mouse':
        genome_index_path = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation.gtf"
    elif genome == 'human':
        genome_index_path = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf"
    
    return genome_index_path,scriptpath_tracks,tracks_dir,macs2_prefix_list,region
    
def visualization_MakeTracks(genome_index_path,scriptpath_tracks,tracks_dir,macs2_prefix_list,region):
    
    peak_list = ""
    for i in range(len(macs2_prefix_list)):
        if i == 0:
            peak_list += macs2_prefix_list[i]
        else:
            peak_list += " "+macs2_prefix_list[i]
    
    jobid = []
    command = "source "+scriptpath_tracks+" "+genome_index_path+" "+peak_list+" "+tracks_dir+" "+region+" "+project_name
    #job = os.popen(command).read().strip().splitlines()
    job = os.popen(command).read().splitlines()
    print(job[1])
    jobid.append(job[1].split(' ')[2])
    
    return jobid

def visualization_PlotTracks(tracks_dir):

    #trackplot = IFrame(tracks_dir+"/"+"genomeTracks.png", width=800, height=800)
    #trackplot = HTML(open(tracks_dir+"/"+"genomeTracks.png").read())
    trackplot = Image(tracks_dir+"/"+"genomeTracks.png")
    
    return trackplot

def visualization_PlotHeatmap(out_peak):

    dirlist = DataFrame(pd.Series(os.listdir(out_peak)))
    dirlist.index = range(len(dirlist))
    print(dirlist)
    
    print("========================================")
    print("Specify index to visualize peak heatmap:(e.g 0)")
    index_files = int(input())    
    
    file = dirlist.loc[index_files,0]
    
    #heatplot = IFrame(out_peak+"/"+file+"/"+file+"_heatmap_TSS.png", width=800, height=800)
    #heatplot = HTML(open(out_peak+"/"+file+"/"+file+"_heatmap_TSS.png").read())
    heatplot = Image(out_peak+"/"+file+"/"+file+"_heatmap_TSS.png")
    
    return heatplot
    
