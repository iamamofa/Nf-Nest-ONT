import argparse
import pandas as pd
import subprocess
import csv
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
args = parser.parse_args()

AAdic= {"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Asx":"B", "Cys":"C", "Glu":"E", "Gln":"Q", "Glx":"Z", "Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K", "Met":"M", "Phe":"F", "Pro":"P", "Ser":"S", "Thr":"T", "Trp":"W", "Tyr":"Y", "Val":"V"}


with open(args.filtered, "r") as f1:
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
                #print(line)
                if count==0:
                    Genelist1+=[word]
                    #print(word)
                if count==1:
                    POSlist1+=[word]
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
                    if word[2:5] not in AAdic or word[5+len(temppos):5+len(temppos)+3] not in AAdic:
                        Genelist1=Genelist1[:-1]
                        POSlist1=POSlist1[:-1]
                        Reflist1=Reflist1[:-1]
                        Altlist1=Altlist1[:-1]
                    #print(word[5+len(temppos):5+len(temppos)+3])
                count+=1
            #print(Genelist1)
            #print(POSlist1)
            #print(Reflist1)
            #print(Altlist1)
            #print(AAlist1)
            #print(AAlist2)
    #print(len(Genelist1))
    #print(len(AAlist2))
    #print(Genelist1)
    #print(AAlist2)
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
    for x in range(len(Genelist1)):
        for y in range(len(depthlist1)):
            if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                #print(depthlist1[y])
                depthlist1up+=[depthlist1[y][2]]
                break
    #print(depthlist1)
    #print(depthlist1up)
    #print(Genelist1)
    #print(POSlist1)
    codondepthlist1=[]
    PfCRTStartlist1=[96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115]
    PfCRTEndlist1=[x+y for x,y in zip([96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115],[90,268,172,132,71,75,82,50,56,92,44,54,76])]
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
    #sadcount=0
    for x in range(len(Genelist1)):
        for y in range(len(depthlist1)):
            if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                #print(depthlist1[y])
                #sadcount+=1
                tempcodondepthlist1=0
                if Genelist1[x] == "PfCRT": 
                    for z in range(len(PfCRTStartlist1)):
                        if PfCRTStartlist1[z]<int(depthlist1[y-1][2])<PfCRTEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y-1][2])
                        if PfCRTStartlist1[z]<int(depthlist1[y][2])<PfCRTEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y][2])
                        if PfCRTStartlist1[z]<int(depthlist1[y+1][2])<PfCRTEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y+1][2])
                        #print(tempcodondepthlist1)
                    codondepthlist1+=[int(tempcodondepthlist1/3)]
                    break
                if Genelist1[x] == "DHFR":
                    for z in range(len(DHFRStartlist1)):
                        #print(type(DHFRStartlist1[z]))
                        #print((depthlist1[y-1][2]))
                        if DHFRStartlist1[z]<int(depthlist1[y-1][2])<DHFREndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y-1][2])
                        if DHFRStartlist1[z]<int(depthlist1[y][2])<DHFREndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y][2])
                        if DHFRStartlist1[z]<int(depthlist1[y+1][2])<DHFREndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y+1][2])
                    codondepthlist1+=[int(tempcodondepthlist1/3)]
                    break
                if Genelist1[x] == "K13": 
                    for z in range(len(K13Startlist1)):
                        if K13Startlist1[z]<int(depthlist1[y-1][2])<K13Endlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y-1][2])
                        if K13Startlist1[z]<int(depthlist1[y][2])<K13Endlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y][2])
                        if K13Startlist1[z]<int(depthlist1[y+1][2])<K13Endlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y+1][2])
                    codondepthlist1+=[int(tempcodondepthlist1/3)]
                    break
                if Genelist1[x] == "PfMDR1":
                    for z in range(len(PfMDR1Startlist1)):
                        if PfMDR1Startlist1[z]<int(depthlist1[y-1][2])<PfMDR1Endlist1[z]:
                                tempcodondepthlist1+=int(depthlist1[y-1][2])
                        if PfMDR1Startlist1[z]<int(depthlist1[y][2])<PfMDR1Endlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y][2])
                        if PfMDR1Startlist1[z]<int(depthlist1[y+1][2])<PfMDR1Endlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y+1][2]) 
                    codondepthlist1+=[int(tempcodondepthlist1/3)]
                    break
                if Genelist1[x] == "DHPS":
                    for z in range(len(PfDHPSStartlist1)):
                        if PfDHPSStartlist1[z]<int(depthlist1[y-1][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y-1][2])
                        if PfDHPSStartlist1[z]<int(depthlist1[y][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y][2])
                        if PfDHPSStartlist1[z]<int(depthlist1[y+1][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y+1][2])
                    codondepthlist1+=[int(tempcodondepthlist1/3)]
                    break 
                if Genelist1[x] == "mitochondrial_genome":
                    for z in range(len(mitochondrialgenomeStartlist1)):
                        if mitochondrialgenomeStartlist1[z]<int(depthlist1[y-1][2])<mitochondrialgenomeEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y-1][2])
                        if mitochondrialgenomeStartlist1[z]<int(depthlist1[y][2])<mitochondrialgenomeEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y][2])
                        if mitochondrialgenomeStartlist1[z]<int(depthlist1[y+1][2])<mitochondrialgenomeEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y+1][2])
                    codondepthlist1+=[int(tempcodondepthlist1/3)]
                    break
                #elif Genelist1[x] == "DHFR": 
                #elif Genelist1[x] == "K13": 
                #elif Genelist1[x] == "PfMDR1": 
                    #break
    #print(codondepthlist1)



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
            #print(count)
            #print(word1)
            #temppos3=""
            #tempgene1=""
            if count==1:
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
                        if word1== Genelist1[x] and word== str(int(POSlist1[x]),) and line.startswith("#")==False:
                            #print(word)
                            #print(word)
                            #filtercount+=1
                            #print(POSlist1[x])
                            #print("test"+POSlist1[x])
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
                                if word.find(wholefilnums[x])!=-1 and [word1,wholefilnumsx[x]] in Geneset1 and [word1,str(int(wholefilnumsx[x]))] not in GenePosSet1:
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

        if QDlist1[x] == "None" and SORlist1[x] == "None" and MQlist1[x] == "None" and MQRankSumlist1[x] == "None":
            Description+="Score is N/A"
        if MQRankSumlist1[x] != "None":
            if MQRankSumlist1[x] >-2:
                Description+=""
            else:
                Description+="MQRankSum is out of range"
        if SORlist1[x] != "None":
            if float(SORlist1[x])<=3:
                Description+=""
            else:
                Description+="SOR is too high "
        #print(QDlist1[filtercount])
        if QDlist1[x] != "None":
            filternumber+=float(QDlist1[x])/80
            if float(QDlist1[x])>=15:
                Description+=""
            else: 
                Description+="QD is too low"
        if MQlist1[x] != "None":
            filternumber+=float(MQlist1[x])/140
            if 35<=float(MQlist1[x]):
                Description+=""
            else:
                Description+="MQ is too low" 
       
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
    print(Genelist2)
    print(POSlist2)
    print(Genelist1)
    print(POSlist1)
    #print(len(Genelist2))
    #print(len(AAPoslist12))
    #print(len(AAlist21))
    #print(len(AAlist22))
    #print(depthlist1up)
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
            wholeCandlist1.append(["PfDHPS",str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
        elif df1Gene[k] == "mitochondrial_genome":
            wholeCandlist1.append([df1Gene[k],str(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
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
            wholeReportlist1.append(["PfDHPS",str(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
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
    d3={"Gene":Genelist2, "BasePOS":POSlist2, "BaseDepth":depthlist2up, "Ref":Reflist2, "Alt":Altlist2, "AAref":AAlist21, "AAalt":AAlist22, "AAPOS":AAPoslist12, "CodonCoverage":codondepthlist2, "AF": AFlist1, "Mutation": mutationlist1, "QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1, "Filter":filterlist1, "FilterDescription":Descriptionlist1, "Candidates":candidateslist1, "Reportable":reportablelist1}
    df3=pd.DataFrame(data=d3)
    #df.to_csv('snpreport4.csv', index=False,) 
    df3.to_csv(args.name+'.csv', index=False,) 
