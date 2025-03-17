import argparse
import pandas as pd
import csv

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-n', dest='name', type=str, help="name of vcf file")

#name='final_SRR6463548_filtered.vcfext.vcf'
args = parser.parse_args()
name=args.name
name2='fixed'+name
name3='fixedPOS'+name
#name3='fixedPOS'+name[:-4]+'.csv'

#print(pd.read_csv(name))

#with open(name, mode='r') as f1:
    #print(f1.read())
#    Lines = f1.readlines() 
#    for line in Lines: 
#        print(len(line.split()), "columns")
#        baseCol=len(line.split())
#        break

#with open(name, mode='r') as f1:
    #print(f1.read())
#    Lines = f1.readlines() 
#    for line in Lines: 
#        if len(line.split()) >  baseCol

with open(name, mode='r') as f1:
    with open(name2, mode='w') as f2:
    #print(f1.read())
        Lines = f1.readlines() 
        for line in Lines: 
            f2.write(line)
            break


with open(name, mode='r') as f1:
    with open(name2, mode='a') as f2:
    #print(f1.read())
        Lines = f1.readlines()
        #testcount=0
        protein=""
        for line in Lines: 
            #testcount+=1
            count=0
            for word in line.split():
                if word.startswith("c."):
                    count+=1
                if word.startswith("n."):
                    count+=1
                if word.startswith("p."):
                    protein=word
            #if protein in line: print("true")
            for x in range((count)):
                count2=-1
                for word in line.split():
                    if word.startswith("c.") == False and word.startswith("n.") == False:
                        f2.write(word+"\t")
                    if word.startswith("c."):
                            count2+=1
                    if word.startswith("n."):
                            count2+=1
                    if x == count2 and protein in line:
                       #print(x)
                        if word.startswith("c."):
                            f2.write(word+"\t"+protein+"\n")
                            break
                        if word.startswith("n."):
                            f2.write(word+"\t"+protein+"\n")
                            break
                    if x == count2 and protein not in line:
                        #print(x)
                        if word.startswith("c."):
                            f2.write(word+"\n")
                            break
                        if word.startswith("n."):
                            f2.write(word+"\n")
                            break
            #if testcount ==2:
            #    break


    with open(name, mode='r') as f1:
        with open(name3, mode='w') as f3:
        #print(f1.read())
            Lines = f1.readlines() 
            for line in Lines: 
                f3.write(line)
                break

    with open(name2, mode='r') as f2:
        with open(name3, mode='a') as f3:
        #print(f1.read())
            Lines = f2.readlines()
            #testcount=0
            start=0
            for line in Lines: 
                #testcount+=1
                start+=1
                tempword=""
                tempcount=0
                count=0
                for word in line.split():
                    tempcount+=1
                    if word.startswith("c.") and "+" not in word and "-" not in word:
                        for character in word: 
                            if character.isdigit():
                                tempword+=character
                    if word.startswith("n.") and "+" not in word and "-" not in word:
                        for character in word: 
                            if character.isdigit():
                                tempword+=character
                    if word.startswith("c.") and "+" in word or "-" in word:
                        start=1
                    if word.startswith("n.") and "+" in word or "-" in word:
                        start=1
                if start>1:
                    for word in line.split():
                        count+=1
                        if count==2:
                            f3.write(tempword+"\t")
                        if count!=2:
                            if count == tempcount:
                                f3.write(word+"\n")
                            else: f3.write(word+"\t")
