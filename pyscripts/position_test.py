import argparse
import pandas as pd
import subprocess
import csv
import sys
def install(name):
    subprocess.call(['pip', 'install', name])

install('xlrd==1.2.0')
import xlrd

AAdic= {"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Asx":"B", "Cys":"C", "Glu":"E", "Gln":"Q", "Glx":"Z", "Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K", "Met":"M", "Phe":"F", "Pro":"P", "Ser":"S", "Thr":"T", "Trp":"W", "Tyr":"Y", "Val":"V"}

with open(sys.argv[1], "r") as f1:
    Genelist1=[]
    POSlist1=[]
    Reflist1=[]
    Altlist1=[]
    AAlist1=[]
    AAlist2=[]
    AAPoslist1=[]
    for line in f1:
        count=0
        if "p." in line:
            for word in line.split():
                #print(word)
                if count==0:
                    Genelist1+=[word]
                if count==1:
                    POSlist1+=[word]
                if count==2:
                    Reflist1+=[word]
                if count==3:
                    Altlist1+=[word]
                if word.startswith("p."):
                    if word[2:5] in AAdic:
                        AAlist1+=[AAdic[word[2:5]]]
                    else:
                        Genelist1=Genelist1[:-1]
                        POSlist1=POSlist1[:-1]
                        Reflist1=Reflist1[:-1]
                        Altlist1=Altlist1[:-1]
                        break
                    if word[-3:] in AAdic:
                        AAlist2+=[AAdic[word[-3:]]]
                        AAPoslist1+=[word[5:-3]]
                        #print(word[5:-3])
                    else:
                        Genelist1=Genelist1[:-1]
                        POSlist1=POSlist1[:-1]
                        Reflist1=Reflist1[:-1]
                        Altlist1=Altlist1[:-1]
                        #AAlist1=AAlist1[:-1]
                        #AAPoslist1=AAPoslist1[:-1]
                        break
                count+=1
print(Genelist1)
print(POSlist1)
            #print(Reflist1)
            #print(Altlist1)
print(AAlist1)
print(AAPoslist1)