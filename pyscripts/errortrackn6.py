from os import listdir
import pandas as pd
import csv
import os
import argparse
import subprocess
def install(name):
    subprocess.call(['pip3', 'install', name])

install('xlrd==1.2.0')
import xlrd

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-n1', dest='nflog', type=str, help="name nextflow nest log")
parser.add_argument('-w1', dest='work', type=str, help="name of work directory")
parser.add_argument('-x1', dest='excel', type=str, help="name of ct file")
args = parser.parse_args()

wholelist1=[["samplename","error", "step", "CTvalue"]]
with open(args.nflog, "r") as logfile:
    samplename=""
    error=""
    step=""
    CTvalue="NA"
    for line in logfile:
        if line.find("error") != -1 and line.find("errortrack") ==-1:
            if line.find("INFO") != -1:
                if line.find("NOTE: Process `"):
                    step=line[line.find("NOTE: Process `")+15:line.find(" (")]
                #print(line[line.find("["):line.find("]")+1])
                #print((line[line.find("]")+1::]))
                #print((line[line.find("]")+1::])[(line[line.find("]")+1::]).find("["):(line[line.find("]")+1::]).find("]")+1])
                workdir=(line[line.find("]")+1::])[(line[line.find("]")+1::]).find("["):(line[line.find("]")+1::]).find("]")+1]
                #print(workdir)
                #print(workdir[1:workdir.find("/")])
                #print(workdir[workdir.find("/")+1:-1])
                workdir1=workdir[1:workdir.find("/")]
                workdir2=workdir[workdir.find("/")+1:-1]
                for x in listdir(args.work+"/"+workdir1):
                    if x.startswith(workdir2):
                        workdir2=x
                with open(args.work+"/"+workdir1+"/"+workdir2+"/.command.sh", "r") as commandsh:
                    for line in commandsh:
                        if line.find("samtools index") !=-1:
                            samplename=line[14::]
                            #print(samplename)
                        if line.find("bbduk.sh") !=-1:
                            samplename=line[line.find("in=")+3:line.find("_R1")]
                        if line.find("AddOrReplaceReadGroups") !=-1:
                            samplename=line[line.find("-I")+3:line.find(".sam")] 
                        if line.find("cat") !=1:
                            samplename=line[line.find(">")+2:line.find("_L001")]
                with open(args.work+"/"+workdir1+"/"+workdir2+"/.command.err", "r") as commanderror:
                    for line in commanderror:
                        #print(line)
                        if line.find("KeyError: 'fields")!=-1:
                            error="noVCF?fields"
                            #print(error)
                            #print(samplename,error)
                        if line.find("KeyError: 'info")!=-1:
                            error="noVCF?info"
                            #print(samplename,error)
                        if line.find("There appear to be different numbers of reads in the paired input files.")!=-1:
                            error="different pair read from trimm"
                            #print(samplename,error)
                        if line.find("AddOrReplaceReadGroups")!=-1:
                            error="java run time error"
                        if line.find("Input is being processed as paired")!=-1:
                            error="fastq error"
                        if line.find("cat: '*R2_001.fastq.gz': No such file or directory")!=-1:
                            error="R2 sequence is missing"
                        if line.find("cat: '*R1_001.fastq.gz': No such file or directory")!=-1:
                            error="R1 sequence is missing"
                #print(samplename)
                #print(error)
                #print(line)
                if os.path.exists(args.excel):
                    if samplename != "":
                        samplename=samplename[:samplename.find("Fxxx0")+1]+"1230"
                    xls=pd.ExcelFile(args.excel)
                    df1 = pd.read_excel(xls, 'TES PCR gel results')
                #print(df1.columns)
                #df1.to_csv("test.csv", index=False)
                #print(df1["Date completed"].tolist())
                #print(df1["Unnamed: 2"].tolist())
                #print(df1['Date completed.1'].tolist())
                #print(df1["Unnamed: 15"].tolist())
                #print(samplename)
                    for x in range(len(df1["Date completed"].tolist())):
                        #print((df1["Date completed"].tolist())[x])
                        #print(samplename)
                        if (df1["Date completed"].tolist())[x]==samplename:
                            #print("Test")
                            CTvalue=(df1["Unnamed: 2"].tolist())[x]
                            wholelist1+=[[samplename,error,step,CTvalue]]
                    for x in range(len(df1["Date completed.1"].tolist())):
                        if (df1["Date completed.1"].tolist())[x]==samplename:
                            #print("Test")
                            #print((df1["Date completed"].tolist())[x])
                            #print(samplename)
                            CTvalue=(df1["Unnamed: 15"].tolist())[x]
                            wholelist1+=[[samplename,error,step,CTvalue]]
                else:
                    wholelist1+=[[samplename,error,step,CTvalue]]
#print(wholelist1)
with open("out.csv", "w", newline="") as f:
    #print(wholelist1)
    writer = csv.writer(f)
    writer.writerows(wholelist1)

                