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
parser.add_argument('-b1', dest='lsam', type=str, help="name of bam file")
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
            #print(Genelist1)
            #print(POSlist1)
            #print(Reflist1)
            #print(Altlist1)
            #print(AAlist1)
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
    for x in range(len(Genelist1)):
        for y in range(len(depthlist1)):
            if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                #print(depthlist1[y])
                depthlist1up+=[depthlist1[y][2]]
                break
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
                if Genelist1[x] == "PfDHPS":
                    for z in range(len(PfDHPSStartlist1)):
                        if PfDHPSStartlist1[z]<int(depthlist1[y-1][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y-1][2])
                        if PfDHPSStartlist1[z]<int(depthlist1[y][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y][2])
                        if PfDHPSStartlist1[z]<int(depthlist1[y+1][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1+=int(depthlist1[y+1][2])
                    codondepthlist1+=[int(tempcodondepthlist1/3)]
                    break 
                if Genelist1[x] == "mitochondrialgenome":
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
    Genelist2=[]
    POSlist2=[]
    AFlist1=[]
    QDlist1=[]
    SORlist1=[]
    MQlist1=[]
    MQRankSumlist1=[]
    filterlist1=[]
    f2= open(args.unfiltered, "r")
    for line in f2:
        count=0
        for word in line.split():
            #print(line)
            #print(Genelist1[x])
            #print(word)
            #if count==0 and word == Genelist1[x]:
                #print(word+"test1111")
            if count==0:
                word1=word
            if count==1:
                filtercount=-1
                for x in range(len(Genelist1)):
                    if word1== Genelist1[x] and word== str(int(POSlist1[x])-1) and line.startswith("#")==False:
                        filtercount+=1
                        Genelist2+=[Genelist1[x]]
                        POSlist2+=[POSlist1[x]]
                        if line.find("AF=")+3 != 2:
                            for k in range(int(line.find("AF="))+3,len(line)):
                                if line[k]==";":
                                    break
                            if line[line.find("AF=")+3:k].find(",") == -1:
                                AFlist1+=[float(line[line.find("AF=")+3:k])]
                            elif line[line.find("AF=")+3:k] .find(",") != -1:
                                AFlist1+=[float(line[line.find("AF=")+3:line.find("AF=")+5])]
                        else:AFlist1+=["None"]
                        if line.find("MQ=")+3 != 2:
                            for k in range(int(line.find("MQ=")+3),len(line)):
                                if line[k]==";":
                                    break
                            MQlist1+=[float(line[line.find("MQ=")+3:k])]
                        else: MQlist1+=["None"]
                        if line.find("MQRankSum=")+10 != 9:
                            for k in range(int(line.find("MQRankSum="))+10,len(line)):
                                if line[k]==";":
                                    break
                            MQRankSumlist1+=[float(line[line.find("MQRankSum=")+10:k])]
                        else:MQRankSumlist1+=["None"]
                        if line.find("SOR=")+4 != 3:
                            for k in range(int(line.find("SOR="))+4,len(line)):
                                if line[k]==";":
                                    break
                            SORlist1+=[float(line[line.find("SOR=")+4:k])]
                        else: SORlist1+=["None"]
                        if line.find("QD=")+3!= 2:
                            for k in range(int(line.find("QD="))+3,len(line)):
                                if line[k]==";":
                                    break
                            QDlist1+=[float(line[line.find("QD=")+3:k])]
                        else:QDlist1+=["None"]
                        #print(filternumber)
            count+=1
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
            filternumber+=float(MQlist1[x])/140
            if 40<float(MQlist1[x])<60:
                Description+="Good,"
            else:
                Description+="Bad,"
        else:Description+="NA,"
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
    for x in range(len(Genelist2)):
        for y in range(len(Genelist1)):
            if Genelist2[x]==Genelist1[y] and POSlist2[x]==POSlist1[y]:
                Reflist2+=[Reflist1[y]]
                Altlist2+=[Altlist1[y]]
                AAlist21+=[AAlist1[y]]
                AAlist22+=[AAlist2[y]]
                depthlist2up+=[depthlist1up[y]]
                AAPoslist12+=[AAPoslist1[y]]
                codondepthlist2+=[codondepthlist1[y]]
                break
    mutationlist1=[]
    for x in range(len(AFlist1)):
        if AFlist1[x]==1:mutationlist1+=["Major"]
        if AFlist1[x]==0.5:mutationlist1+=["Minor"]
        if AFlist1[x]==0:mutationlist1+=["Wildtype"]
        if AFlist1[x]=="None":mutationlist1+=["Wildtype"]

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
    #print(MQlist1)
    #print(MQRankSumlist1)
    #print(len(Genelist1))
    #print(len(Genelist2))
    #print(len(POSlist2))
    #print(len(depthlist2up))
    wholelenlist=[len(POSlist2),len(depthlist2up),len(Reflist2),len(Altlist2),len(AAlist21),len(AAlist22),len(AAPoslist12),len(codondepthlist2),len(AFlist1),len(mutationlist1),len(QDlist1),len(SORlist1),len(MQlist1),len(MQRankSumlist1),len(filterlist1)]
    for x in range(len(wholelenlist)):
        if wholelenlist[x] != len(Genelist2):
            print("order", x, wholelenlist[x], len(Genelist2))
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
    for x in range(len(Genelist2)):
        #print("test")
        if Genelist2[x]=="PfCRT":
            for k in range(len(df1Gene)):
                if df1Gene[k]=="Pfcrt":
                    if AAPoslist12[x]==df1AAPos[k]:
                        if AAlist21[x]==df1RefAA[k]:
                            if AAlist22==df1AltAA[k]:
                                candidateslist1+["candidates"]
                                break
                            else:
                                candidateslist1+=["no"]
                                break
                        else:
                            candidateslist1+=["no"]
                            break
                    else:
                        candidateslist1+=["no"]
                        break
                else:
                    candidateslist1+=["no"]
                    break
        if Genelist2[x]=="PfCRT":
            for k in range(len(df2Gene)):
                if df2Gene[k]=="Pfcrt":
                    #print("test")
                    if AAPoslist12[x]==df2AAPos[k]:
                        if AAlist21[x]==df2RefAA[k]:
                            if AAlist22==df2AltAA[k]:
                                reportablelist1+["reportable"]
                                break
                            else:
                                reportablelist1+=["no"]
                                break
                        else:
                            reportablelist1+=["no"]
                            break
                    else: 
                        reportablelist1+=["no"]
                        break
        if Genelist2[x]=="DHFR":
            #print("test")
            for k in range(len(df1Gene)):
                if df1Gene[k]=="Pfdhfr":
                    #print("test")
                    if AAPoslist12[x]==df1AAPos[k]:
                        if AAlist21[x]==df1RefAA[k]:
                            if AAlist22==df1AltAA[k]:
                                candidateslist1+["candidates"]
                                break
                            else:
                                candidateslist1+=["no"]
                                break
                        else:
                            candidateslist1+=["no"]
                            break
                    else:
                        candidateslist1+=["no"]
                        break
                else:
                    candidateslist1+=["no"]
                    break
        if Genelist2[x]=="DHFR":
            for k in range(len(df2Gene)):
                if df2Gene[k]=="Pfdhfr":
                    if AAPoslist12[x]==df2AAPos[k]:
                        if AAlist21[x]==df2RefAA[k]:
                            if AAlist22==df2AltAA[k]:
                                reportablelist1+["reportable"]
                                break
                            else:
                                reportablelist1+=["no"]
                                break
                        else:
                            reportablelist1+=["no"]
                            break
                    else: 
                        reportablelist1+=["no"]
                        break
        if Genelist2[x]=="K13":
            for k in range(len(df1Gene)):
                if df1Gene[k]=="PfK13":
                    if AAPoslist12[x]==df1AAPos[k]:
                        if AAlist21[x]==df1RefAA[k]:
                            if AAlist22==df1AltAA[k]:
                                candidateslist1+["candidates"]
                                break
                            else:
                                candidateslist1+=["no"]
                                break
                        else:
                            candidateslist1+=["no"]
                            break
                    else:
                        candidateslist1+=["no"]
                        break
                else:
                    candidateslist1+=["no"]
                    break
        if Genelist2[x]=="K13":
            for k in range(len(df2Gene)):
                if df2Gene[k]=="Pfk13":
                    if AAPoslist12[x]==df2AAPos[k]:
                        if AAlist21[x]==df2RefAA[k]:
                            if AAlist22==df2AltAA[k]:
                                reportablelist1+["reportable"]
                                break
                            else:
                                reportablelist1+=["no"]
                                break
                        else:
                            reportablelist1+=["no"]
                            break
                    else: 
                        reportablelist1+=["no"]
                        break
        if Genelist2[x]=="PfMDR1":
            for k in range(len(df1Gene)):
                if df1Gene[k]=="Pfmdr1":
                    if AAPoslist12[x]==df1AAPos[k]:
                        if AAlist21[x]==df1RefAA[k]:
                            if AAlist22==df1AltAA[k]:
                                candidateslist1+["candidates"]
                                break
                            else:
                                candidateslist1+=["no"]
                                break
                        else:
                            candidateslist1+=["no"]
                            break
                    else:
                        candidateslist1+=["no"]
                        break
                else:
                    candidateslist1+=["no"]
                    break
        if Genelist2[x]=="PfMDR1":
            for k in range(len(df2Gene)):
                if df2Gene[k]=="Pfmdr1":
                    if AAPoslist12[x]==df2AAPos[k]:
                        if AAlist21[x]==df2RefAA[k]:
                            if AAlist22==df2AltAA[k]:
                                reportablelist1+["reportable"]
                                break
                            else:
                                reportablelist1+=["no"]
                                break
                        else:
                            reportablelist1+=["no"]
                            break
                    else: 
                        reportablelist1+=["no"]
                        break
        if Genelist2[x]=="PfDHPS":
            for k in range(len(df1Gene)):
                if df1Gene[k]=="Pfdhps":
                    if AAPoslist12[x]==df1AAPos[k]:
                        if AAlist21[x]==df1RefAA[k]:
                            if AAlist22==df1AltAA[k]:
                                candidateslist1+["candidates"]
                                break
                            else:
                                candidateslist1+=["no"]
                                break
                        else:
                            candidateslist1+=["no"]
                            break
                    else:
                        candidateslist1+=["no"]
                        break
                else:
                    candidateslist1+=["no"]
                    break
        if Genelist2[x]=="PfDHPS":
            for k in range(len(df2Gene)):
                if df2Gene[k]=="Pfdhps":
                    if AAPoslist12[x]==df2AAPos[k]:
                        if AAlist21[x]==df2RefAA[k]:
                            if AAlist22==df2AltAA[k]:
                                reportablelist1+["reportable"]
                                break
                            else:
                                reportablelist1+=["no"]
                                break
                        else:
                            reportablelist1+=["no"]
                            break
                    else: 
                        reportablelist1+=["no"]
                        break
        if Genelist2[x]=="mitochondrialgenome":
            for k in range(len(df1Gene)):
                if df1Gene[k]=="mitochondrialgenome":
                    if AAPoslist12[x]==df1AAPos[k]:
                        if AAlist21[x]==df1RefAA[k]:
                            if AAlist22==df1AltAA[k]:
                                candidateslist1+["candidates"]
                                break
                            else:
                                candidateslist1+=["no"]
                                break
                        else:
                            candidateslist1+=["no"]
                            break
                    else:
                        candidateslist1+=["no"]
                        break
                else:
                    candidateslist1+=["no"]
                    break
        if Genelist2[x]=="mitochondrialgenome":
            for k in range(len(df2Gene)):
                if df2Gene[k]=="mitochondrialgenome":
                    if AAPoslist12[x]==df2AAPos[k]:
                        if AAlist21[x]==df2RefAA[k]:
                            if AAlist22==df2AltAA[k]:
                                reportablelist1+["reportable"]
                                break
                            else:
                                reportablelist1+=["no"]
                                break
                        else:
                            reportablelist1+=["no"]
                            break
                    else: 
                        reportablelist1+=["no"]
                        break
    #print(Genelist2)
    #print(reportablelist1)     
    #print(candidateslist1)           
    d3={"Gene":Genelist2, "BasePOS":POSlist2, "BaseDepth":depthlist2up, "Ref":Reflist2, "Alt":Altlist2, "AAref":AAlist21, "AAalt":AAlist22, "AAPOS":AAPoslist12, "CodonCoverage":codondepthlist2, "AF": AFlist1, "Mutation": mutationlist1, "QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1, "Filter":filterlist1, "FilterDescription":Descriptionlist1, "Candidates":candidateslist1, "Reportable":reportablelist1}
    df3=pd.DataFrame(data=d3)
    #df.to_csv('snpreport4.csv', index=False,) 
    df3.to_csv(args.name+'.csv', index=False,) 
