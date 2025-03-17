import argparse
import pandas as pd
import subprocess
import csv
import itertools
def install(name):
    subprocess.call(['pip', 'install', name])

install('xlrd==1.2.0')
import xlrd
parser = argparse.ArgumentParser(description='name')
parser.add_argument('-v1', dest='unfiltered', type=str, help="name of unfilterd merged vcf file")
parser.add_argument('-v2', dest='filtered', type=str, help="name of filtered merged vcf file")
parser.add_argument('-b1', dest='bam', type=str, help="name of bam file")
parser.add_argument('-o1', dest='name', type=str, help="name of output")
parser.add_argument('-e1', dest='candidates', type=str, help="name of candidates")
parser.add_argument('-e2', dest='voi', type=str, help="name of variants of interest")
parser.add_argument('-f1', dest='fasta', type=str, help="name of fasta file")
args = parser.parse_args()

AAdic= {"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Asx":"B", "Cys":"C", "Glu":"E", "Gln":"Q", "Glx":"Z", "Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K", "Met":"M", "Phe":"F", "Pro":"P", "Ser":"S", "Thr":"T", "Trp":"W", "Tyr":"Y", "Val":"V", "STOP":"STOP"}

currentcodon=""
keywords = [''.join(i) for i in itertools.product(["T","G","A","C"], repeat = 3)]
codondic={}
for x in keywords:
    currentcodon=""
    for y in x:
        currentcodon+=y
    #print(currentcodon)
    if currentcodon == "TTT" or currentcodon == "TTC":
        #print("True")
        codondic[currentcodon]="Phe"
    if currentcodon == "TTA" or currentcodon == "TTG": 
        codondic[currentcodon]="Leu"
    if currentcodon == "CTT" or currentcodon =="CTC" or currentcodon =="CTA" or currentcodon =="CTG":
        codondic[currentcodon]="Leu"
    if currentcodon == "ATT" or currentcodon =="ATC" or currentcodon =="ATA":
        codondic[currentcodon]="Ile"
    if currentcodon == "ATG":
        codondic[currentcodon]="Met"
    if currentcodon == "GTT" or currentcodon =="GTC" or currentcodon =="GTA" or currentcodon =="GTG":
        codondic[currentcodon]="Val"
    if currentcodon == "TCT" or currentcodon =="TCC" or currentcodon =="TCA" or currentcodon =="TCG":
        codondic[currentcodon]="Ser"
    if currentcodon == "CCT" or currentcodon =="CCC" or currentcodon =="CCA" or currentcodon =="CCG": 
        codondic[currentcodon]="Pro"
    if currentcodon == "ACT" or currentcodon =="ACC" or currentcodon =="ACA" or currentcodon =="ACG":
        codondic[currentcodon]="Thr"
    if currentcodon == "GCT" or currentcodon =="GCC" or currentcodon =="GCA" or currentcodon =="GCG":
        codondic[currentcodon]="Ala"
    if currentcodon == "TAT" or currentcodon =="TAC":
        codondic[currentcodon]="Tyr"
    if currentcodon == "CAT" or currentcodon =="CAC":
        codondic[currentcodon]="His"
    if currentcodon == "CAA" or currentcodon =="CAG":
        codondic[currentcodon]="Gln"
    if currentcodon == "AAT" or currentcodon =="AAC":
        codondic[currentcodon]="Asn"
    if currentcodon == "AAA" or currentcodon =="AAG":
        codondic[currentcodon]="Lys"
    if currentcodon == "GAT" or currentcodon =="GAC":
        codondic[currentcodon]="Asp"
    if currentcodon == "GAA" or currentcodon =="GAG":
        codondic[currentcodon]="Glu"
    if currentcodon == "TGT" or currentcodon== "TGC":
        #print("I am correct")
        codondic[currentcodon]="Cys"
    if currentcodon == "TGG":
        codondic[currentcodon]="Trp"
    if currentcodon == "AGA" or currentcodon =="AGG" or currentcodon =="CGT" or currentcodon =="CGC" or currentcodon =="CGA" or currentcodon =="CGG":
        codondic[currentcodon]="Arg"
    if currentcodon == "AGT" or currentcodon == "AGC":
        codondic[currentcodon]="Ser"
    if currentcodon == "GGT" or currentcodon =="GGC" or currentcodon =="GGA" or currentcodon =="GGG":
        codondic[currentcodon]="Gly"
    if currentcodon == "TGA" or currentcodon == "TAA" or currentcodon == "TAG":
        codondic[currentcodon]="STOP"
    #print(currentcodon)
    #print(currentcodon)
    #if currentcodon=="TGT" or currentcodon=="TGC":
    #    print("I am correct")
    #print(codondic[currentcodon])

#print(codondic)


#PfCRTStartlist1=[96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115]
#PfCRTEndlist1=[x+y for x,y in zip([96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115],[90,268,172,132,71,75,82,50,56,92,44,54,76])]
#DHFRStartlist1=[1]
#DHFREndlist1=[x+y for x,y in zip([1],[1827])]
#K13Startlist1=[1]
#K13Endlist1=[x+y for x,y in zip([1],[2181])]
#PfMDR1Startlist1=[1]
#PfMDR1Endlist1=[x+y for x,y in zip([1],[4260])]
#PfDHPSStartlist1=[1,313,2302]
#PfDHPSEndlist1=[x+y for x,y in zip([1,313,2302],[135,1868,115])]
#mitochondrialgenomeStartlist1=[734,1933,3492]
#mitochondrialgenomeEndlist1=[x+y for x,y in zip([734,1933,3492],[839,1538,1130])]

with open(args.filtered, "r") as f1:
    #with open(args.fasta, "r") as f2:
    Genelist1=[]
    POSlist1=[]
    Reflist1=[]
    Altlist1=[]
    AAlist1=[]
    AAlist2=[]
    AAPoslist1=[]
    currentgene1=""
    currentpos1=""
    currentAminoacid=""
    currentposlen1=0
    for line in f1:
        if "p." in line:
            count=0
            for word in line.split():
                #print(line)
                if count==0:
                    Genelist1+=[word]
                    currentgene1=word
                    #print(word)
                if count==1:
                    POSlist1+=[word]
                    #currentposlen1=len(POSlist1)
                    #if currentgene1 == "PfCRT": 
                    #    for z in range(len(PfCRTStartlist1)):
                    #        if PfCRTStartlist1[z]<int(currentpos1)<PfCRTEndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "DHFR":
                    #    for z in range(len(DHFRStartlist1)):
                    #        #print(type(DHFRStartlist1[z]))
                            #print((depthlist1[y-1][2]))
                    #        if DHFRStartlist1[z]<int(currentpos1)<DHFREndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "K13": 
                    #    for z in range(len(K13Startlist1)):
                    #        if K13Startlist1[z]<int(currentpos1)<K13Endlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "PfMDR1":
                    #    for z in range(len(PfMDR1Startlist1)):
                    #        if PfMDR1Startlist1[z]<int(currentpos1)<PfMDR1Endlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "PfDHPS":
                    #    for z in range(len(PfDHPSStartlist1)):
                    #        if PfDHPSStartlist1[z]<int(currentpos1)<PfDHPSEndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "mitochondrial_genome":
                    #    for z in range(len(mitochondrialgenomeStartlist1)):
                    #        if mitochondrialgenomeStartlist1[z]<int(currentpos1)<mitochondrialgenomeEndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if len(POSlist1)==currentposlen1:
                        #print("test")
                    #   Genelist1=Genelist1[:-1]
                    #    break
                    #print(word)
                if count==2:
                    Reflist1+=[word]
                if count==3:
                    Altlist1+=[word]
                if word.startswith("p."):
                    #print(word)
                    if word[2:5] in AAdic:
                        AAlist1+=[AAdic[word[2:5]]]
                    temppos=""
                    for x in word:
                        #print(word)
                        if x.isdigit():
                            #print(x)
                            temppos+=x
                    #print(res)
                    if word[-1:] == "*":
                        AAlist2+=["*"]
                        AAPoslist1+=[temppos]
                        break
                    if word[5+len(temppos):5+len(temppos)+3] in AAdic:
                        AAlist2+=[AAdic[word[5+len(temppos):5+len(temppos)+3]]]
                        AAPoslist1+=[temppos]
                    if (word[2:5] not in AAdic and word[-1:] == "*") or (word[2:5] not in AAdic and word[5+len(temppos):5+len(temppos)+3] in AAdic):
                        AAlist2=AAlist2[:-1]
                        AAPoslist1=AAPoslist1[:-1]
                    if word[2:5] not in AAdic or word[5+len(temppos):5+len(temppos)+3] not in AAdic:
                        Genelist1=Genelist1[:-1]
                        POSlist1=POSlist1[:-1]
                        Reflist1=Reflist1[:-1]
                        Altlist1=Altlist1[:-1]
                    #print(word[5+len(temppos):5+len(temppos)+3])
                count+=1
        elif ("p." not in line and "c." in line) or ("p." not in line and "n." in line):
            count=0
            currentAminoacid=""
            for word in line.split():
                #print(line)
                if count==0:
                    Genelist1+=[word]
                    currentgene1=word
                    #print(word)
                if count==1:
                    POSlist1+=[word]
                    currentpos1=word
                    #print(currentpos1)
                if count==2:
                    Reflist1+=[word]
                if count==3:
                    Altlist1+=[word]
                count+=1
            linecount=0
            #print(currentgene1)
            #print(currentpos1)
            with open(args.fasta, "r") as f2:
                for line in f2:
                    #print(currentgene1)
                    #print(linecount)
                    if linecount==1:
                        linecount+=1
                    if line.find(currentgene1)!=-1:
                        #print(line)
                        #print(currentpos1)
                        linecount+=1
                    if linecount==2:
                        #print(currentgene1)
                        #print(currentpos1)
                        currentAminoacid = line[int(currentpos1)-2:int(currentpos1)+1]
            #print(currentAminoacid)
            #print(currentAminoacid)
            if currentAminoacid in codondic:
                #print("TRUE")
                if codondic[currentAminoacid] in AAdic:
                    #print("TRUE")
                    #print(AAdic[codondic[line[int(currentpos1)-2:int(currentpos1)+1]]])
                    AAlist1+=[AAdic[codondic[currentAminoacid]]]
                    AAlist2+=[AAdic[codondic[currentAminoacid]]]
                    AAPoslist1+=[int(currentpos1)/3]
                    #line[int(currentpos1)-2:int(currentpos1)+1]
    #print(len(Genelist1))
    #print(len(POSlist1))
    #print(len(Reflist1))
    #print(len(Altlist1))
    #print(len(AAlist1))
    #print(len(AAPoslist1))
    #print(len(AAlist2))
    #print(len(Genelist1))
    #print(len(AAlist2))
    #print(len(Genelist1))
    #print(Genelist1)
    #print(AAlist2)
    #print(POSlist1)
    process=subprocess.Popen(['samtools', 'depth', '-a', args.name+ "_SR.bam"], stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    reader = csv.reader(stdout.decode('ascii').splitlines(),
                            delimiter='\t', skipinitialspace=True)
    depthlist1=[]
    for row in reader:
        depthlist1+=[row]
    depthlist1up=[]
    #for x in range(len(Genelist1)):
    #    if Genelist1[x] in depthlist1:
    #        print(Genelist1[x])
    #    if POSlist1[x] in depthlist1:
    #        print(POSlist1[x])
    #print(depthlist1)
    #print(depthlist1)
    #print(Genelist1)
    #print(POSlist1)
    #print(Genelist1)
    depthpair1=[]
    for x in range(len(depthlist1)):
        depthpair1+=[[depthlist1[x][0],depthlist1[x][1]]]

    for x in range(len(Genelist1)):
        if [Genelist1[x], POSlist1[x]] not in depthpair1:
            depthlist1up+=[None]

        else:
            for y in range(len(depthlist1)):
                if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                    #print(depthlist1[y])
                    #print(depthlist1[y][2])
                    depthlist1up+=[depthlist1[y][2]]
                    break
    #print(len(depthlist1up))
    #print(depthpair1)
    #print(depthlist1)
    #print(len(depthlist1up))
    #print(depthlist1up)
    #print(Genelist1)
    #print(POSlist1)
    codondepthlist1=[]
    #sadcount=0
    for x in range(len(Genelist1)):
        if [Genelist1[x],POSlist1[x]] not in depthpair1:
                #print("True")
                codondepthlist1+=[None]
                #break
        else:
            for y in range(len(depthlist1)):
                if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                    #print(Genelist1[x],POSlist1[x])
                    #print(depthlist1[y])
                    #sadcount+=1
                    tempcodondepthlist1=0
                    if Genelist1[x] == "PfCRT": 
                        #for z in range(len(PfCRTStartlist1)):
                            #if PfCRTStartlist1[z]<int(depthlist1[y-1][2])<PfCRTEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if PfCRTStartlist1[z]<int(depthlist1[y][2])<PfCRTEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if PfCRTStartlist1[z]<int(depthlist1[y+1][2])<PfCRTEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                            #print(tempcodondepthlist1)
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "DHFR":
                        #for z in range(len(DHFRStartlist1)):
                            #print(type(DHFRStartlist1[z]))
                            #print((depthlist1[y-1][2]))
                            #if DHFRStartlist1[z]<int(depthlist1[y-1][2])<DHFREndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if DHFRStartlist1[z]<int(depthlist1[y][2])<DHFREndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if DHFRStartlist1[z]<int(depthlist1[y+1][2])<DHFREndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "K13": 
                        #for z in range(len(K13Startlist1)):
                            #if K13Startlist1[z]<int(depthlist1[y-1][2])<K13Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if K13Startlist1[z]<int(depthlist1[y][2])<K13Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if K13Startlist1[z]<int(depthlist1[y+1][2])<K13Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "PfMDR1":
                        #for z in range(len(PfMDR1Startlist1)):
                            #if PfMDR1Startlist1[z]<int(depthlist1[y-1][2])<PfMDR1Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if PfMDR1Startlist1[z]<int(depthlist1[y][2])<PfMDR1Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if PfMDR1Startlist1[z]<int(depthlist1[y+1][2])<PfMDR1Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2]) 
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                    if Genelist1[x] == "DHPS":
                        #for z in range(len(PfDHPSStartlist1)):
                            #if PfDHPSStartlist1[z]<int(depthlist1[y-1][2])<PfDHPSEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if PfDHPSStartlist1[z]<int(depthlist1[y][2])<PfDHPSEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if PfDHPSStartlist1[z]<int(depthlist1[y+1][2])<PfDHPSEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "mitochondrial_genome":
                        #for z in range(len(mitochondrialgenomeStartlist1)):
                            #if mitochondrialgenomeStartlist1[z]<int(depthlist1[y-1][2])<mitochondrialgenomeEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if mitochondrialgenomeStartlist1[z]<int(depthlist1[y][2])<mitochondrialgenomeEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if mitochondrialgenomeStartlist1[z]<int(depthlist1[y+1][2])<mitochondrialgenomeEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                #elif Genelist1[x] == "DHFR": 
                #elif Genelist1[x] == "K13": 
                #elif Genelist1[x] == "PfMDR1": 
                    #break
    #print(len(Genelist1))
    #print(len(codondepthlist1))



#PfCRT   1   3399    PfCRT   .   +   96  3191    0   13  90,268,172,132,71,75,82,50,56,92,44,54,76,  96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115,
#MT  1   5967    COXIII  .   -   734 1573    0   1   839,    734,
#MT  1   5967    COL .   +   1933    3471    0   1   1538,   1933,
#MT  1   5967    CYTOb   .   +   3492    4622    0   1   1130,   3492,
#PfDHFR  1   1827    PfDHFR  .   +   1   1827    0   1   1827,   1,
#PfDHPS  1   2417    PfDHPS  .   +   1   2417    0   3   135,1868,115,   1,313,2302,
#PfK13   1   2181    PfK13   .   +   1   2181    0   1   2181,   1,
#PfMDR1  1   4260    PfMDR1  .   +   1   4260    0   1   4260,   1,

    
    #print(len(Genelist1))
    #print(len(POSlist1))
    #print(len(depthlist1up))
    #print(depthlist1up)
    #with open(args.unfiltered, "r") as f2:
    #    for line in f2:
    #        print(line)
   # print(len(Genelist1))
    #print(POSlist1)
    #print(POSlist1)
    #print(Genelist1)
    Genelist2=[]
    POSlist2=[]
    AFlist1=[]
    QDlist1=[]
    SORlist1=[]
    MQlist1=[]
    MQRankSumlist1=[]
    filterlist1=[]
    f2= open(args.unfiltered, "r")
    f3= open(args.filtered, "r")
    wholefilnumsx=[]
    wholefilnums=[]
    for line in f3: 
        count=0
        tempnum=""
        #print(line)
        if line.find("p.")!=-1 and line.startswith("#")==False:
            for word in line.split():
                count+=1
                if count==2:
                    tempnum=word
                #print(tempnum)
                    #print(word)
                #print(tempnum)
                if word.startswith("c."):
                    if word.find("*")!=-1:
                        wholefilnumsx+=[tempnum]
                        #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                        wholefilnums+=[word[word.find("c.*"):word.find(tempnum)+len(tempnum)]]
                    elif word.find("*")==-1:
                        #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                        wholefilnumsx+=[tempnum]
                        wholefilnums+=[word[word.find("c."):word.find(tempnum)+len(tempnum)]]
                if word.startswith("p."):
                    #print(word[5+len(temppos):5+len(temppos)+3])
                    #print(word[5+len(temppos):5+len(temppos)+3])
                    temppos=""
                    for x in word:
                        #print(word)
                        if x.isdigit():
                            #print(x)
                            temppos+=x
                    if word[2:5] not in AAdic or word[5+len(temppos):5+len(temppos)+3] not in AAdic:
                        #print("working")
                        if word[2:5] in AAdic and word[-1:] == "*":
                            break
                        wholefilnumsx=wholefilnumsx[:-1]
                        wholefilnums=wholefilnums[:-1]
    #print(wholefilnumsx)
            #print(tempnum)
    #print(wholefilnums)
    #print(len(Genelist1))
    #print(len(wholefilnums))
    #print(wholefilnums)
    #print(wholefilnumsx)
    #print(wholefilnums)
    #print(wholefilnums)
    #print(wholefilnumsx)
    #print(Genelist1)
    #print(POSlist1)
    Geneset1=[]
    for k in range(len(Genelist1)):
        Geneset1+=[[Genelist1[k],POSlist1[k]]]
    GenePosSet1=[]
    for line in f2:
        count=0
        for word in line.split():
            #print(line)
            #print(Genelist1[x])
            #print(word)
            #if count==0 and word == Genelist1[x]:
                #print(word+"test1111")
            if count==0:
                word1=""
                if word.startswith("#")==False:
                    word1=word
                    #print(word1)
            #print(word1)
            #print(count)
            #print(word1)
            #temppos3=""
            #tempgene1=""
            if count==1:
                #print(word1)
                #print(word1)
                if word1!="":
                    #filtercount=-1
                    #print(word1)
                    #print(word)
                    #print(POSlist1)
                    #print(Genelist1)
                    #print(word1)
                    for x in range(len(Genelist1)):
                        #print(str(int(POSlist1[x])-1))
                        #if word1== Genelist1[x]:print(word1,word)
                        #print(POSlist1)
                        #print(word)
                        #if word.startswith("c."):
                            #temppos2=""
                            #for x in word:
                                #print(word)
                                #if x.isdigit():
                                    #print(x)
                                    #temppos2+=x
                            #print(temppos2)
                        #print(word1)
                        #print(Genelist1[x])
                        #print(POSlist1[x])
                        #print(word)
                        #print(Genelist1[x],POSlist1[x])
                        #print(word1,word)
                        #print(POSlist1)
                        #print(type(POSlist1[x]))
                        #if word1 == Genelist1[x]:print("no")
                        #print(word)
                        #print(POSlist1[x])
                        #if word == POSlist1[x]: print("yes")
                        if word1== Genelist1[x] and word== POSlist1[x] and line.startswith("#")==False:
                            #print(word)
                            #print(word)
                            #print(word)
                            #filtercount+=1
                            #print(POSlist1[x])
                            #print("test"+POSlist1[x])
                            #print(Genelist1[x],POSlist1[x])
                            Genelist2+=[Genelist1[x]]
                            POSlist2+=[POSlist1[x]]
                            GenePosSet1+=[[Genelist1[x],POSlist1[x]]]
                            #temppos3=POSlist1[x]
                            if line.find(";AF=") != -1:
                                #print(line[int(line.find("AF="))+3:len(line)].find(";"))
                                #print(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")])
                                #print(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")].isdigit())
                                #tempAF=""
                                countAF=0
                                for k in (line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")]):
                                    if k.isdigit() == False and k !=".":
                                        AFlist1+=[line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")]]
                                        break
                                    #print(len(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")]))
                                    #print("COUNTAF"+str(countAF))
                                    if countAF == len(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")])-1:
                                        AFlist1+=[line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")]]
                                        break
                                    countAF+=1
                                #else:AFlist1+=[line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")]]
                                #for k in range(int(line.find("AF="))+3,line[int(line.find("AF="))+3:len(line)].find(";")):
                                #    if line[k]==";":
                                #        break
                                #    if line[line.find("AF=")+3:k].find(",") == -1:
                                #        AFlist1+=[float(line[line.find("AF=")+3:k])]
                                #    elif line[line.find("AF=")+3:k].find(",") != -1:
                                #        AFlist1+=[float(line[line.find("AF=")+3:line.find("AF=")+5])]
                            else: AFlist1+=["None"]
                            #print(AFlist1)
                            if line.find("MQ=") != -1:
                            #if line.find("MQ=")+3 != 2:
                            #    for k in range(int(line.find("MQ=")+3),len(line)):
                            #        if line[k]==";":
                            #            break
                            #        MQlist1+=[float(line[line.find("MQ=")+3:k])]
                                MQlist1+=[line[int(line.find("MQ="))+3:int(line.find("MQ="))+3+line[int(line.find("MQ="))+3:len(line)].find(";")]]
                            else: MQlist1+=["None"]
                            if line.find("MQRankSum=") != -1:
                            #    for k in range(int(line.find("MQRankSum="))+10,len(line)):
                            #        if line[k]==";":
                            #           break
                            #        MQRankSumlist1+=[float(line[line.find("MQRankSum=")+10:k])]
                                #print(int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10])
                                #print(int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";"))
                                #print(line[int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";")])
                                MQRankSumlist1+=[line[int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";")]]
                            else:MQRankSumlist1+=["None"]
                            if line.find("SOR=") != -1:
                                #for k in range(int(line.find("SOR="))+4,len(line)):
                                #    if line[k]==";":
                                #        break
                                #    SORlist1+=[float(line[line.find("SOR=")+4:k])]
                                SORlist1+=[line[int(line.find("SOR="))+4:int(line.find("SOR="))+4+line[int(line.find("SOR="))+4:len(line)].find(";")]]
                            else: SORlist1+=["None"]
                            if line.find("QD=") != -1:
                                #for k in range(int(line.find("QD="))+3,len(line)):
                                #    if line[k]==";":
                                #        break
                                #    QDlist1+=[float(line[line.find("QD=")+3:k])]
                                QDlist1+=[line[int(line.find("QD="))+3:int(line.find("QD="))+3+line[int(line.find("QD="))+3:len(line)].find(";")]]
                            else: QDlist1+=["None"]
                            #print(filternumber)
                            break
            #print(len(QDlist1))
            #print(GenePosSet1)
            #print(AFlist1)
            if count==7:
                if word1 != "":
                    #print(word)
                    if line.startswith("#")==False:
                        if word.startswith("AC=") or word.startswith("AB="):
                            #print(word[word.find("c."):word.find("c.")+10])
                            #print(line)
                            for x in range(len(wholefilnums)):
                                #print(x)
                                #print(wholefilnumsx[count])
                                #print(POSlist2)
                                #print(wholefilnumsx[count2])
                                #print([word1,str(int(wholefilnumsx[x]))])
                                #print(word)
                                #find part of c.138 in c.1381 in different gene.....
                                #print(line.find(wholefilnums[x]))
                                #if word.find(wholefilnums[x]) != -1:
                                    #print(wholefilnums[x])
                                    #print(word.find(wholefilnums[x]))
                                    #print(len(wholefilnums[x]))
                                    #print(wholefilnums[x])
                                    #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])])
                                    #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])-len(wholefilnumsx[x])-1])
                                if word.find(wholefilnums[x])!=-1 and str(word.find(wholefilnums[x])+len(wholefilnums[x])+1).isdigit()==False and [word1,wholefilnumsx[x]] in Geneset1 and [word1,str(int(wholefilnumsx[x]))] not in GenePosSet1:
                                    #print(word1)
                                    #print(wholefilnumsx[x])
                                    #print(wholefilnums[x])
                                    #print(wholefilnumsx[x])
                                    #print(word1)
                                    #print(word1)
                                    #print(wholefilnumsx[x])
                                    Genelist2+=[word1]
                                    #print(wholefilnumsx[x])
                                    POSlist2+=[wholefilnumsx[x]]
                                    #"QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1,
                                    AFlist1+=[AFlist1[-1]]
                                    QDlist1+=[QDlist1[-1]]
                                    SORlist1+=[SORlist1[-1]]
                                    MQlist1+=[MQlist1[-1]]
                                    MQRankSumlist1+=[MQRankSumlist1[-1]]
                                    GenePosSet1+=[[word1,str(int(wholefilnumsx[x]))]]
                                #print(Genelist1[count])
                                #print(wholefilnumsx[count])
            count+=1
    #print(len(wholefilnums))
    #print(len(Genelist2))
    #print(len(QDlist1))
    #print(len(Genelist2))
    #print(POSlist1)
    #print(POSlist2)
    #print(Genelist1)
    #print(Genelist2)
    Descriptionlist1=[]
    for x in range(len(QDlist1)):
        filternumber=0
        Description=""
        #print(QDlist1[filtercount])
        if QDlist1[x] != "None":
            filternumber+=float(QDlist1[x])/80
            if float(SORlist1[x])>30:
                Description+="Good,"
            else:
                Description+="Bad,"
        else:Description+="NA,"
        if SORlist1[x] != "None":
            if float(SORlist1[x])<3:
                Description+="Good,"
            else:
                Description+="Bad,"
        else:Description+="NA,"
        if MQlist1[x] != "None":
            if float(MQlist1[x])<70:
                filternumber+=float(MQlist1[x])/140
            else:filternumber+=1/2
            if 40<float(MQlist1[x])<60:
                Description+="Good,"
            else:
                Description+="Bad,"
        else:Description+="NA,"
        #print(MQRankSumlist1)
        if MQRankSumlist1[x] != "None":
            if -5<abs(float(MQRankSumlist1[x]))<5:
                Description+="Good"
            else:
                Description+="Bad"
        else:Description+="NA"
        Descriptionlist1+=[Description]
        #print(Description)
        filterlist1+=[filternumber]
    #print(filterlist1)
                        #print(word+"test2222")
    #for x in range(len(Genelist1)):
        #print(Genelist1[x])
    #    for line in f2:
            #count=0
            #print(line)
    #        for word in line.split():
                #if count==0:
                    #print(word + "test")
    #            if Genelist1[x] == word:
    #                print(Genelist1[x] +"test222")
                    #    print(line)
                    #    break
            #    count+=1
                #if count==1:
                #    if POSlist1[x] == word:
                #        print(line)

    #print(len(Genelist1))
    #print(len(depthlist1up))
    #(sadcount)
    #print(len(AAPoslist1))
    #print(len(codondepthlist1))
    Reflist2=[]
    Altlist2=[]
    AAlist21=[]
    AAlist22=[]
    depthlist2up=[]
    AAPoslist12=[]
    codondepthlist2=[]
    #print(Genelist2)
    #print(POSlist2)
    #print(Genelist1)
    #print(POSlist1)
    #print(len(Genelist2))
    #print(len(POSlist2))
    #print(len(AAPoslist12))
    #print(len(AAlist21))
    #print(len(AAlist22))
    #print(len(depthlist1up))
    #print(len(Genelist1))
    #print(len(AAlist2))
    #print(len(Reflist1))
    count=0
    for x in range(len(Genelist2)):
        for y in range(len(Genelist1)):
            if  Genelist2[x]==Genelist1[y] and POSlist2[x]==POSlist1[y]:
                count+=1
                Reflist2+=[Reflist1[y]]
                Altlist2+=[Altlist1[y]]
                AAlist21+=[AAlist1[y]]
                AAlist22+=[AAlist2[y]]
                depthlist2up+=[depthlist1up[y]]
                AAPoslist12+=[AAPoslist1[y]]
                codondepthlist2+=[codondepthlist1[y]]
                break
    #Geneset2=[]
    #for k in range(len(Genelist1)):
    #    Geneset2+=[[Genelist1[k],POSlist1[k]]]
    #for x in range(len(Genelist2)):
    #    for y in range(len(Genelist1)):
    #        if [Genelist2[x],POSlist2[x]] not in Geneset2:
    #            print(Genelist2[x],POSlist2[x])
    #for y in POSlist2:
    #    if y not in POSlist1:
    #        print(y)
    #print(len(Genelist2))
    #print(count)
    mutationlist1=[]
    #print(AFlist1)
    for x in range(len(AFlist1)):
        if AFlist1[x]!="None":
            if float(AFlist1[x])==1.0:mutationlist1+=["Major"]
            if float(AFlist1[x])==0.5:mutationlist1+=["Minor"]
            if float(AFlist1[x])==0.0:mutationlist1+=["Wildtype"]
        if AFlist1[x]=="None":mutationlist1+=["Wildtype"]

    #print(len(wholefilnumsx))
    #print(len(wholefilnums))
    #print(len(Genelist1))
    #print(len(Genelist2))
    #print(len(POSlist2))
    #print(len(AAPoslist12))
    #print(len(AAlist21))
    #print(len(AAlist22))

    #print(len(AFlist1))
    #print(len(QDlist1))
    #print(len(MQlist1))
    #print(len(SORlist1))
    #print(len(MQRankSumlist1))
    #print(len(filterlist1))
    #print(len(AAlist22))
    #print(len(AAPoslist12))
    #d1={"Gene":Genelist1, "POS":POSlist1, "Ref":Reflist1, "Alt":Altlist1, "AAref":AAlist1, "AAalt":AAlist2, "Depth":depthlist1up}
    #df=pd.DataFrame(data=d1)
    #d2={"Gene":Genelist2, "BasePOS":POSlist2, "Ref":Reflist2, "Alt":Altlist2, "AAref":AAlist21, "AAalt":AAlist22, "Depth":depthlist2up, "AF": AFlist1, "Mutation": mutationlist1, "QD": QDlist1}
    #df2=pd.DataFrame(data=d2)
    #print(Genelist2)
    #print(POSlist2)
    #print(Genelist1)
    #print(POSlist1)
    #print(MQlist1)
    #print(MQRankSumlist1)
    #print(len(Genelist1))
    #print(len(Genelist2))
    #print(len(POSlist2))
    #print(len(depthlist2up))
    #wholelenlist=[len(POSlist2),len(depthlist2up),len(Reflist2),len(Altlist2),len(AAlist21),len(AAlist22),len(AAPoslist12),len(codondepthlist2),len(AFlist1),len(mutationlist1),len(QDlist1),len(SORlist1),len(MQlist1),len(MQRankSumlist1),len(filterlist1)]
    #for x in range(len(wholelenlist)):
    #    if wholelenlist[x] != len(Genelist2):
    #        print("order", x, wholelenlist[x], len(Genelist2))
    xls=pd.ExcelFile(args.candidates)
    df1=pd.read_excel(xls)
    df1Gene=df1['Gene'].tolist()
    df1RefAA=df1['RefAA'].tolist()
    df1AltAA=df1['AltAA'].tolist()
    df1AAPos=df1['AAPos'].tolist()
    #df1.to_csv("test1.csv", index=False)
    df2 = pd.read_csv(args.voi)
    df2Gene=df2['Gene'].tolist()
    df2RefAA=df2['RefAA'].tolist()
    df2AltAA=df2['AltAA'].tolist()
    df2AAPos=df2['AAPos'].tolist()
    #print(df2RefAA)
    #print(df2Gene)
    #print(df2RefAA)
    candidateslist1=[]
    reportablelist1=[]
    #print(df1Gene)
    #print(df1AAPos)
    #print(df1RefAA)
    #print(df1AltAA)
    wholeAAlist1=[]
    wholeCandlist1=[]
    wholeReportlist1=[]
    #print(len(Genelist2))
    #print(len(Genelist1))
    #print(len(AAPoslist12))
    #print(len(AAlist21))
    #print(len(AAlist22))
    #print(Genelist2)
    for x in range(len(Genelist2)):
        #print(len(Genelist2))
        #print(len(Genelist1))
        #print(len(AAPoslist12))
        #print(len(AAlist21))
        #print(len(AAlist22))
        wholeAAlist1.append([Genelist2[x],AAPoslist12[x],AAlist21[x],AAlist22[x]])
    for k in range(len(df1Gene)):
        if df1Gene[k] == "Pfcrt":
            wholeCandlist1.append(["PfCRT",str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
        elif df1Gene[k] == "Pfdhfr":
            wholeCandlist1.append(["DHFR",str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
        elif df1Gene[k] == "PfK13":
            wholeCandlist1.append(["K13",str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
        elif df1Gene[k] == "Pfmdr1":
            wholeCandlist1.append(["PfMDR1",str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
        elif df1Gene[k] == "Pfdhps":
            wholeCandlist1.append(["DHPS",str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
        elif df1Gene[k] == "MT":
            wholeCandlist1.append(["mitochondrial_genome",str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
    for k in range(len(df2Gene)):
        if df2Gene[k] == "Pfcrt":
            wholeReportlist1.append(["PfCRT",str(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
        elif df2Gene[k] == "Pfdhfr":
            wholeReportlist1.append(["DHFR",str(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
        elif df2Gene[k] == "PfK13":
            wholeReportlist1.append(["K13",str(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
        elif df2Gene[k] == "Pfmdr1":
            wholeReportlist1.append(["PfMDR1",str(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
        elif df2Gene[k] == "Pfdhps":
            wholeReportlist1.append(["DHPS",str(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
        elif df2Gene[k] == "mitochondrial_genome":
            wholeReportlist1.append([df2Gene[k],str(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
    #print(wholeAAlist1)
    #print(wholeReportlist1)
    #print(['PfCRT', 74, 'M', 'I'])
    #candidateslist1=[]
    for x in ((wholeAAlist1)):
        #print(x)
        if x in wholeCandlist1:
            candidateslist1+=["candidates"]
        else: candidateslist1+=["no"]
    #    print(x)
        if x in wholeReportlist1:
            reportablelist1+=["reportable"]
        else:reportablelist1+=["no"]
        #print("test")
        #print(Genelist2)
        #print(Genelist2[x])
        #print(AAPoslist12[x])
        #print(AAlist21[x])
        #print(AAlist22[x])
                            #print(tempAlt1[tempcount2])
                            #tempcount2+=1
                    #f3.write(templine+"\n")
            #if tempcount!=len(GeneNamelist1):
            #    f3.write(">"+tempGene1+"\n")

    #print(codondic)
    #print(tempdic)
    filteredGenelist2=[]
    filteredPoslist2=[]
    filtereddepthlist2up=[]
    filteredReflist2=[]
    filteredAltlist2=[]
    filteredAAlist21=[]
    filteredAAlist22=[]
    filteredAAPoslist12=[]
    filteredcodondepthlist2=[]
    filteredAFlist1=[]
    filteredmutationlist1=[]
    filteredQDlist1=[]
    filteredSORlist1=[]
    filteredMQlist1=[]
    filteredMQRankSumlist1=[]
    filteredfilterlist1=[]
    filteredDescriptionlist1=[]
    filteredcandidateslist1=[]
    filteredreportablelist1=[]
    #print(df1AAPos)
    #print(df2AAPos)
    #print(Genelist1)
    for x in range(len(AAPoslist12)):
        if float(AAPoslist12[x]) in df1AAPos or float(AAPoslist12[x]) in df2AAPos:
            #print(AAPoslist12[x])
            filteredGenelist2+=[Genelist2[x]]
            filteredPoslist2+=[POSlist2[x]]
            filtereddepthlist2up+=[depthlist2up[x]]
            filteredReflist2+=[Reflist2[x]]
            filteredAltlist2+=[Altlist2[x]]
            filteredAAlist21+=[AAlist21[x]]
            filteredAAlist22+=[AAlist22[x]]
            filteredAAPoslist12+=[AAPoslist12[x]]
            filteredcodondepthlist2+=[codondepthlist2[x]]
            filteredAFlist1+=[AFlist1[x]]
            filteredmutationlist1+=[mutationlist1[x]]
            filteredQDlist1+=[QDlist1[x]]
            filteredSORlist1+=[SORlist1[x]]
            filteredMQlist1+=[MQlist1[x]]
            filteredMQRankSumlist1+=[MQRankSumlist1[x]]
            filteredfilterlist1+=[filterlist1[x]]
            filteredDescriptionlist1+=[Descriptionlist1[x]]
            filteredcandidateslist1+=[candidateslist1[x]]
            filteredreportablelist1+=[reportablelist1[x]]

    #print(Genelist1)
    #print(Genelist2)
    #print(reportablelist1)     
    #print(candidateslist1)
    #print(len(Genelist2))           
    #print(len(POSlist2))
    #print(len(depthlist2up))
    #print(len(Reflist2))
    #print(len(Altlist2))
    #print(len(AAlist21))
    #print(len(AAlist22))
    #print(len(AAPoslist12))
    #print(len(codondepthlist2))
    #print(len(AFlist1))
    #print(len(mutationlist1))
    #print(len(QDlist1))
    #print(len(SORlist1))
    #print(len(MQlist1))
    #print(len(MQRankSumlist1))
    #print(len(filterlist1))
    #print(len(Descriptionlist1))
    #print(len(candidateslist1))
    #print(len(reportablelist1))
    #print(len())
    #print(len())
    d3={"Gene":filteredGenelist2, "BasePOS":filteredPoslist2, "BaseDepth":filtereddepthlist2up, "Ref":filteredReflist2, "Alt":filteredAltlist2, "AAref":filteredAAlist21, "AAalt":filteredAAlist22, "AAPOS":filteredAAPoslist12, "CodonCoverage":filteredcodondepthlist2, "AF": filteredAFlist1, "Mutation": filteredmutationlist1, "QD": filteredQDlist1, "SOR":filteredSORlist1, "MQ":filteredMQlist1, "MQRankSum":filteredMQRankSumlist1, "Filter":filteredfilterlist1, "FilterDescription":filteredDescriptionlist1, "Candidates":filteredcandidateslist1, "Reportable":filteredreportablelist1}
    df3=pd.DataFrame(data=d3)
    #df.to_csv('snpreport4.csv', index=False,) 
    df3.to_csv(args.name+'.csv', index=False,) 


            #print(templine)


            #print(tempPOS1)

                        #    print(l)
            #for x in line:
            #if linecount==1:
            #    linecount+=1
            #if line.find(currentgene1)!=-1:
                #print(line)
                #print(currentpos1)
            #    linecount+=1
            #if linecount==2:
                #print(currentgene1)
                #print(currentpos1)
            #    currentAminoacid = line[int(currentpos1)-2:int(currentpos1)+1]
    #print(f3)
    #tempdic={}


    PfCRTStartlist1=[96,366,812,1158,1443,1638,1812,2020,2208,2413,2699,2891,3115]
    PfCRTEndlist1=[x+y for x,y in zip([96,366,812,1158,1443,1638,1812,2020,2208,2413,2699,2891,3115],[89,266,170,131,71,75,80,50,56,92,44,54,76])]
    DHFRStartlist1=[1]
    DHFREndlist1=[x+y for x,y in zip([1],[1827])]
    K13Startlist1=[1]
    K13Endlist1=[x+y for x,y in zip([1],[2181])]
    PfMDR1Startlist1=[1]
    PfMDR1Endlist1=[x+y for x,y in zip([1],[4260])]
    PfDHPSStartlist1=[1,313,2302]
    PfDHPSEndlist1=[x+y for x,y in zip([1,313,2302],[135,1868,115])]
    mitochondrialgenomeStartlist1=[734,1933,3492]
    mitochondrialgenomeEndlist1=[x+y for x,y in zip([734,1933,3492],[839,1538,1130])]
    AAtranslatedic={}
    with open(args.fasta, "r") as f4:
        tempGene1=""
        for line in f4:
            if count!=1:
                count=0
            if count==1:
                count+=1
            if line[1::].strip()=="PfCRT":
                tempGene1="PfCRT"
                count+=1
            if count==2 and tempGene1=="PfCRT":
                for x in range(len(PfCRTStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(PfCRTStartlist1[x]-1,PfCRTEndlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="DHFR":
                tempGene1="DHFR"
                count+=1
            if count==2 and tempGene1=="DHFR":
                for x in range(len(DHFRStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(DHFRStartlist1[x]-1,DHFREndlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="K13":
                tempGene1="K13"
                count+=1
            if count==2 and tempGene1=="K13":
                for x in range(len(K13Startlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(K13Startlist1[x]+1,K13Endlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="PfMDR1":
                tempGene1="PfMDR1"
                count+=1
            if count==2 and tempGene1=="PfMDR1":
                for x in range(len(PfMDR1Startlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(PfMDR1Startlist1[x]-1,PfMDR1Endlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="DHPS":
                tempGene1="DHPS"
                count+=1
            if count==2 and tempGene1=="DHPS":
                for x in range(len(PfDHPSStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(PfDHPSStartlist1[x]-1,PfDHPSEndlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="mitochondrial_genome":
                tempGene1="mitochondrial_genome"
                count+=1
            if count==2 and tempGene1=="mitochondrial_genome":
                for x in range(len(mitochondrialgenomeStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(mitochondrialgenomeStartlist1[x]-2,mitochondrialgenomeEndlist1[x]-2):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
    #print(AAtranslatedic)
    #print(AAtranslatedic)
        #for x in range(len(DHFRStartlist1)):
        #for x in range(len(K13Startlist1)):
        #for x in range(len(PfMDR1Startlist1)):
        #for x in range(len(PfDHPSStartlist1)):
        #for x in range(len(mitochondrialgenomeStartlist1)):

    GeneNamelist1=["PfCRT","DHFR","DHPS","PfMDR1","mitochondrial_genome", "K13"]
    with open(args.fasta, "r") as f2:
        wholetranAA2=[]
        with open("translatedfasta.fa", "w") as f3:
            count=0
            #tempGene1=""
            #tempGene2=""
            #wholelinelist1=[]
            tempGene2=""
            for line in f2:
                tempcount=0
                #if tempGene1!="":
                #    tempGene2=tempGene1
                tempGene1=""
                templine=""
                for k in GeneNamelist1:
                    if line.find(k)!=-1:
                        tempPOS1=[]
                        tempAlt1=[]
                        tempGene1=line[1::].strip()
                        #print(tempGene1)
                        for l in range(len(filteredGenelist2)):
                            #print((l))
                            #print((tempGene1))
                            #print(l == tempGene1)
                            if filteredGenelist2[l] == tempGene1:
                                tempPOS1+=[int(filteredPoslist2[l])]
                                tempAlt1+=[filteredAltlist2[l]]
                for k in GeneNamelist1:
                    if line.find(k)==-1:
                        tempcount+=1
                #if tempGene1!="":
                #    print(tempGene1)
                #print(tempGene2)
                if tempPOS1!=[]:
                    templine=""
                    tempcount2=0
                    if tempcount==len(GeneNamelist1):
                        for x in range(len(line.strip())):
                            if (x+1) not in tempPOS1:
                                templine+=line.strip()[x]
                            else:
                                templine+=tempAlt1[tempcount2]
                                if (tempGene2,x+1) in AAtranslatedic:
                                    wholetranAA2+=[AAtranslatedic[tempGene2,x+1]]
                                #print(tempAlt1[tempcount2])
                                tempcount2+=1
                        f3.write(templine+"\n")
                if tempPOS1==[] and tempcount==len(GeneNamelist1):
                    for x in range(len(line.strip())):
                            if (x+1) not in tempPOS1:
                                templine+=line.strip()[x]
                    f3.write(templine+"\n")
                if tempcount!=len(GeneNamelist1):
                    f3.write(">"+tempGene1+"\n")
                    tempGene2=tempGene1
    #print(wholetranAA2)

    AAtranslatedic2={}
    with open("translatedfasta.fa", "r") as f4:
        tempGene1=""
        for line in f4:
            if count!=1:
                count=0
            if count==1:
                count+=1
            if line[1::].strip()=="PfCRT":
                tempGene1="PfCRT"
                count+=1
            if count==2 and tempGene1=="PfCRT":
                for x in range(len(PfCRTStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(PfCRTStartlist1[x]-1,PfCRTEndlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="DHFR":
                tempGene1="DHFR"
                count+=1
            if count==2 and tempGene1=="DHFR":
                for x in range(len(DHFRStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(DHFRStartlist1[x]-1,DHFREndlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="K13":
                tempGene1="K13"
                count+=1
            if count==2 and tempGene1=="K13":
                for x in range(len(K13Startlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(K13Startlist1[x]+1,K13Endlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="PfMDR1":
                tempGene1="PfMDR1"
                count+=1
            if count==2 and tempGene1=="PfMDR1":
                for x in range(len(PfMDR1Startlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(PfMDR1Startlist1[x]-1,PfMDR1Endlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="DHPS":
                tempGene1="DHPS"
                count+=1
            if count==2 and tempGene1=="DHPS":
                for x in range(len(PfDHPSStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(PfDHPSStartlist1[x]-1,PfDHPSEndlist1[x]+1):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1
            if line[1::].strip()=="mitochondrial_genome":
                tempGene1="mitochondrial_genome"
                count+=1
            if count==2 and tempGene1=="mitochondrial_genome":
                for x in range(len(mitochondrialgenomeStartlist1)):
                    tempcount=0
                    tempcodon=""
                    for y in range(mitochondrialgenomeStartlist1[x]-2,mitochondrialgenomeEndlist1[x]-2):
                        tempcount+=1
                        #print(tempGene1,y)
                        if tempcount==1:
                            tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                            AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                        if tempcount==2:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==3:
                            AAtranslatedic2[tempGene1,y+1]=tempcodon
                        if tempcount==4:
                            if len(line.strip()[y::])>2:
                                tempcodon=AAdic[codondic[line.strip()[y:y+3]]]
                                AAtranslatedic2[tempGene1,y+1]=AAdic[codondic[line.strip()[y:y+3]]]
                                tempcount=1

    with open("translatedfasta.fa", "r") as f4:
        #print("test")
        count=0
        wholetranGene1=[]
        wholetranPOS1=[]
        wholetranAA1=[]
        for line in f4:
            tempcount=0
            #if tempGene1!="":
            #    tempGene2=tempGene1
            tempGene1=""
            templine=""
            for k in GeneNamelist1:
                    if line.find(k)!=-1:
                        tempPOS1=[]
                        tempAlt1=[]
                        tempdic={}
                        tempGenelist1=[]
                        tempGene1=line[1::].strip()
                        #print(tempGene1)
                        for l in range(len(filteredGenelist2)):
                            #print((l))
                            #print((tempGene1))
                            #print(l == tempGene1)
                            if filteredGenelist2[l] == tempGene1:
                                #print(POSlist2[l])
                                tempdic[filteredPoslist2[l]]=filteredGenelist2[l]
                                tempPOS1+=[int(filteredPoslist2[l])]
                                tempAlt1+=[filteredAltlist2[l]]
            for k in GeneNamelist1:
                if line.find(k)==-1:
                    tempcount+=1
            #if tempGene1!="":
            #    print(tempGene1)
            #print(tempGene2)
            #print(temPOS1)
            if tempPOS1!=[] and tempdic!={}:
                templine=""
                tempcount2=0
                if tempcount==len(GeneNamelist1):
                    for x in range(len(line.strip())):
                        #print(x)
                        if (x+1) in tempPOS1:
                            #print(line.strip()[x:x+3])
                            #print(codondic[line.strip()[x:x+3]])
                            wholetranGene1+=[tempdic[str(x+1)]]
                            wholetranPOS1+=[x+1]
                            wholetranAA1+=[AAtranslatedic2[tempdic[str(x+1)],x+1]]
                            #print(tempdic[str(x+1)])
                            #print(x+1)
                            #print(AAdic[codondic[line.strip()[x:x+3]]])

    d3={"Gene":wholetranGene1, "BasePOS":wholetranPOS1, "previousAA": wholetranAA2, "customAA":wholetranAA1}
    df3=pd.DataFrame(data=d3)
    df3.to_csv('translate1.csv', index=False,) 
                            


with open(args.filtered, "r") as f1:
#with open(args.fasta, "r") as f2:
    tempGene1=""
    tempPOS1=""
    Genelist1=[]
    POSlist1=[]
    GenePOSSet1=[]
    for line in f1:
        if line.find("ANN[*].EFFECT")==-1:
            count=0
            for word in line.split():
                #print(line)
                if count==0:
                    tempGene1=word
                    #print(word)
                if count==1:
                    tempPOS1=word
                if count==2:
                    if [tempGene1,tempPOS1] not in GenePOSSet1:
                        Genelist1+=[tempGene1]
                        POSlist1+=[tempPOS1]
                        GenePOSSet1+=[[tempGene1,tempPOS1]]
                count+=1



    filterlist1=[]
    f2= open(args.unfiltered, "r")
    f3= open(args.filtered, "r")
    wholefilnumsx=[]
    wholefilnums=[]
    wholenumset=[]
    for line in f3: 
        count=0
        tempnum=""
        tempGene2=""
        #print(line)
        if line.find("ANN[*].EFFECT")==-1:
            for word in line.split():
                count+=1
                if count==1:
                    tempGene2=word
                if count==2:
                    tempnum=word
                #print(tempnum)
                    #print(word)
                #print(tempnum)
                if word.startswith("c.") or word.startswith("n."):
                    #print(word)
                    if word.find("*")!=-1 and [tempGene2,tempnum] not in wholenumset:
                        wholefilnumsx+=[tempnum]
                        #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                        wholefilnums+=[word[word.find("c.*"):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]
                    if word.find("*")==-1 and [tempGene2,tempnum] not in wholenumset:
                        #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                        wholefilnumsx+=[tempnum]
                        wholefilnums+=[word[word.find("c."):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]
                    if word.find("*")!=-1 and [tempGene2,tempnum] not in wholenumset:
                        wholefilnumsx+=[tempnum]
                        wholefilnums+=[word[word.find("n.*"):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]
                    if word.find("*")==-1 and [tempGene2,tempnum] not in wholenumset:
                        wholefilnumsx+=[tempnum]
                        wholefilnums+=[word[word.find("n."):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]

#print(wholefilnumsx)
#print(wholefilnums)


Genelist2=[]
POSlist2=[]
GenePosSet1=[]
Geneset1=[]
customAA1=[]
for k in range(len(Genelist1)):
    Geneset1+=[[Genelist1[k],POSlist1[k]]]

f2= open(args.unfiltered, "r")
for line in f2:
    tempcount=0
    count=0
    for word in line.split():
        #print(line)
        #print(Genelist1[x])
        #print(word)
        #if count==0 and word == Genelist1[x]:
            #print(word+"test1111")
        if count==0:
            word1=""
            if word.startswith("#")==False:
                word1=word
                #print(word1)
        #print(word1)
        #print(count)
        #print(word1)
        #temppos3=""
        #tempgene1=""
        if count==1:
            #print(word1)
            #print(word1)
            if word1!="":
                #filtercount=-1
                #print(word1)
                #print(word)
                #print(POSlist1)
                #print(Genelist1)
                for x in range(len(Genelist1)):
                    #print(str(int(POSlist1[x])-1))
                    #if word1== Genelist1[x]:print(word1,word)
                    #print(POSlist1)
                    #print(word)
                    #if word.startswith("c."):
                        #temppos2=""
                        #for x in word:
                            #print(word)
                            #if x.isdigit():
                                #print(x)
                                #temppos2+=x
                        #print(temppos2)
                    #print(word1)
                    #print(Genelist1[x])
                    #print(POSlist1[x])
                    #print(word)
                    #print(Genelist1[x],POSlist1[x])
                    #print(word1,word)
                    #print(POSlist1)
                    #print(type(POSlist1[x]))
                    #if word1 == Genelist1[x]:print("no")
                    #print(word)
                    #print(POSlist1[x])
                    #if word == POSlist1[x]: print("yes")
                    if word1== Genelist1[x] and word== POSlist1[x] and line.startswith("#")==False:
                        #print(word1)
                        #print(word)
                        #print(word)
                        #print(word)
                        #filtercount+=1
                        #print(POSlist1[x])
                        #print("test"+POSlist1[x])
                        #print(Genelist1[x],POSlist1[x])
                        Genelist2+=[Genelist1[x]]
                        POSlist2+=[POSlist1[x]]
                        GenePosSet1+=[[Genelist1[x],POSlist1[x]]]
                        #print(line.find("CUSTOM&annotation")+18)
                        #print(line)
                        if line.find("CUSTOM&annotation")!=-1:
                            #print(line.find("CUSTOM&annotation"))
                            #print(line[2::])
                            #print(line[line.find("CUSTOM&annotation")+18::])
                            #print(line[line.find("CUSTOM&annotation")+18::].find("|"))
                            customAA1+=[line[line.find("CUSTOM&annotation")+18:line.find("CUSTOM&annotation")+18+line[line.find("CUSTOM&annotation")+18::].find("|")]]
                        else:
                            customAA1+=["None"]
                        #temppos3=POSlist1[x]
        #print(len(QDlist1))
        #print(GenePosSet1)
        #print(AFlist1)
        if count==7:
            if word1 != "":
                #print(word)
                if line.startswith("#")==False:
                    if word.startswith("AC=") or word.startswith("AB="):
                        #print(word[word.find("c."):word.find("c.")+10])
                        #print(line)
                        for x in range(len(wholefilnums)):
                            #print(x)
                            #print(wholefilnumsx[count])
                            #print(POSlist2)
                            #print(wholefilnumsx[count2])
                            #print([word1,str(int(wholefilnumsx[x]))])
                            #print(word)
                            #find part of c.138 in c.1381 in different gene.....
                            #print(line.find(wholefilnums[x]))
                            #if word.find(wholefilnums[x]) != -1:
                                #print(wholefilnums[x])
                                #print(word.find(wholefilnums[x]))
                                #print(len(wholefilnums[x]))
                                #print(wholefilnums[x])
                                #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])])
                                #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])-len(wholefilnumsx[x])-1])
                            if word.find(wholefilnums[x])!=-1 and str(word.find(wholefilnums[x])+len(wholefilnums[x])+1).isdigit()==False and [word1,wholefilnumsx[x]] in Geneset1 and [word1,str(int(wholefilnumsx[x]))] not in GenePosSet1:
                                #print(word1)
                                #print(wholefilnumsx[x])
                                #print(wholefilnums[x])
                                #print(wholefilnumsx[x])
                                #print(word1)
                                #print(word1)
                                #print(wholefilnumsx[x])
                                Genelist2+=[word1]
                                #print(wholefilnumsx[x])
                                POSlist2+=[wholefilnumsx[x]]
                                #"QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1,
                                GenePosSet1+=[[word1,str(int(wholefilnumsx[x]))]]
                                if line.find("CUSTOM&annotation")!=-1:
                            #print(line.find("CUSTOM&annotation"))
                            #print(line[2::])
                            #print(line[line.find("CUSTOM&annotation")+18::])
                            #print(line[line.find("CUSTOM&annotation")+18::].find("|"))
                                    customAA1+=[line[line.find("CUSTOM&annotation")+18:line.find("CUSTOM&annotation")+18+line[line.find("CUSTOM&annotation")+18::].find("|")]]
                                else:
                                    customAA1+=["None"]
                            #print(Genelist1[count])
                            #print(wholefilnumsx[count])
        count+=1

#print(Genelist1)
#print(POSlist1)
#print(Genelist2)
#print(POSlist2)
#print(customAA1)
d3={"Gene":Genelist2, "BasePOS":POSlist2, "customAA":customAA1}
df3=pd.DataFrame(data=d3)
df3.to_csv('custom1.csv', index=False,) 
