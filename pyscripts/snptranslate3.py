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
parser.add_argument('-l1', dest='lowcov', type=str, help="name of lowcov", default=None)
args = parser.parse_args() 

dict1={}
with open(args.lowcov, "r") as l1:
    lowcov1 = csv.reader(l1, delimiter=' ')
    for line in lowcov1:
        tempword1=""
        tempnum1=""
        count=0
        for word in line[0]:
            if word=="," or word=="\n":
                count+=1
            #print(word)
            #print(count)
            if count==0:
                tempword1+=word
            if count==1:
                tempnum1+=word
        #print(tempword1)
        #print(tempnum1[1::])
        dict1[tempword1]=tempnum1[1::]

#if file doesn't exist  os.mkdir(args.directory2)
if os.path.isdir(args.directory2) == False:
    os.mkdir(args.directory2)

for filename in os.listdir(args.directory):
    if "DS_Store" not in filename:
        count=0
        with open(os.path.join(args.directory, filename), "r") as f1:
            for line in f1:
                list1=line.split(",")
                if list1[0]!="Gene" and list1[2]!="NA":
                    #print(filename)
                    #print(line)
                    if int(dict1[list1[0]])>int(list1[2]):
                        if count==0:
                            with open(os.path.join(args.directory2, filename[:-4]+"2.csv"), "w") as f1:
                                f1.write("Gene,BasePOS,BaseDepth,Ref,Alt,AAref,AAalt,AAPOS,CodonCoverage,AF,Mutation,QD,SOR MQ,MQRankSum,Filter,FilterDescription,Candidates,Reportables\n")
                                for item in list1:
                                    if item != list1[-1]:
                                        f1.write(item+",")
                                    if item == list1[-1]:
                                        f1.write(item)
                        if count>0:
                             with open(os.path.join(args.directory2, filename[:-4]+"2.csv"), "a") as f1:
                                for item in list1:
                                    if item != list1[-1]:
                                        f1.write(item+",")
                                    if item == list1[-1]:
                                        f1.write(item)
                        count+=1
                if list1[0]!="Gene" and list1[2]=="NA":
                    if count==0:
                        with open(os.path.join(args.directory2, filename[:-4]+"2.csv"), "w") as f1:
                            f1.write("Gene,BasePOS,BaseDepth,Ref,Alt,AAref,AAalt,AAPOS,CodonCoverage,AF,Mutation,QD,SOR MQ,MQRankSum,Filter,FilterDescription,Candidates,Reportables\n")
                            for item in list1:
                                if item != list1[-1]:
                                    f1.write(item+",")
                                if item == list1[-1]:
                                    f1.write(item)
                    if count>0:
                         with open(os.path.join(args.directory2, filename[:-4]+"2.csv"), "a") as f1:
                            for item in list1:
                                if item != list1[-1]:
                                    f1.write(item+",")
                                if item == list1[-1]:
                                    f1.write(item)
                    count+=1
                #print(list1[0])
                #print(list1[2])

