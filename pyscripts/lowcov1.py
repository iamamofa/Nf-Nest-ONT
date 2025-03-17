import argparse
import pandas as pd
import subprocess
import csv
from collections import OrderedDict
import json
import os

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-d1', dest='directory', type=str, help="name of directory", default=None)
parser.add_argument('-d2', dest='directory2', type=str, help="name of directory", default=None)
args = parser.parse_args() 

#print("19USxxxxx2HB3PfDxxx0_S95_L001_R2_001.fastq.gz" in os.listdir(args.directory2))

contlist1=[]
for filename in os.listdir(args.directory2):
    contlist1+=[filename[:-16]]


dict1={}
for filename in os.listdir(args.directory):
    #print(filename[:-3]+"fastq.gz")
    if filename[:-4] in contlist1:
        #print(filename)
        with open(os.path.join(args.directory, filename), "r") as f1:
            for line in f1:
                #print(line)
                if line.startswith("DHPS"):
                    count=0
                    tempword=""
                    for word in line:
                        if word==",":
                            count+=1
                        if count==2 and word!=",":
                            tempword+=word
                    if "DHPS" in dict1:
                        if tempword<dict1["DHPS"]:
                            dict1["DHPS"]=tempword
                    else:dict1["DHPS"]=tempword
                if line.startswith("DHFR"):
                    count=0
                    tempword=""
                    for word in line:
                        if word==",":
                            count+=1
                        if count==2 and word!=",":
                            tempword+=word
                    if "DHFR" in dict1:
                        if tempword<dict1["DHFR"]:
                            dict1["DHFR"]=tempword
                    else:dict1["DHFR"]=tempword
                if line.startswith("PfCRT"):
                    count=0
                    tempword=""
                    for word in line:
                        if word==",":
                            count+=1
                        if count==2 and word!=",":
                            tempword+=word
                    if "PfCRT" in dict1:
                        if tempword<dict1["PfCRT"]:
                            dict1["PfCRT"]=tempword
                    else:dict1["PfCRT"]=tempword
                if line.startswith("PfMDR1"):
                    count=0
                    tempword=""
                    for word in line:
                        if word==",":
                            count+=1
                        if count==2 and word!=",":
                            tempword+=word
                    if "PfMDR1" in dict1:
                        if tempword<dict1["PfMDR1"]:
                            dict1["PfMDR1"]=tempword
                    else:dict1["PfMDR1"]=tempword
                if line.startswith("K13"):
                    count=0
                    tempword=""
                    for word in line:
                        if word==",":
                            count+=1
                        if count==2 and word!=",":
                            tempword+=word
                    if "K13" in dict1:
                        if tempword<dict1["K13"]:
                            dict1["K13"]=tempword
                    else:dict1["K13"]=tempword
                if line.startswith("mitochondrial_genome"):
                    count=0
                    tempword=""
                    for word in line:
                        if word==",":
                            count+=1
                        if count==2 and word!=",":
                            tempword+=word
                    if "mitochondrial_genome" in dict1:
                        if tempword<dict1["mitochondrial_genome"]:
                            dict1["mitochondrial_genome"]=tempword
                    else:dict1["mitochondrial_genome"]=tempword


with open("lowcov1.csv", "w") as l1:
    for item in dict1:
        l1.write(item+","+dict1[item]+"\n")


