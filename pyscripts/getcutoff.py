import pandas as pd
import os
import argparse
import glob


#os.chdir("/Users/subinpark/Downloads/nextflow_cloud/test_nextflow_NeST_docker/local/angola_coverage_cut1/samtoolcoverage/")
parser = argparse.ArgumentParser(description='name')
parser.add_argument('-c', dest='coveragefile',type=str, help="name of coverage files")
args = parser.parse_args()


df=pd.DataFrame()
for file in os.listdir("."):
    if file.endswith("coverageall.txt"):
        f=open(file)
        df=df.append(pd.DataFrame(pd.read_table(f)))
#df=df[(df['coverage'] != 0) & (df.numreads)]
df.groupby("#rname")[['numreads','coverage']].describe()
q1=df.groupby("#rname")[['coverage']].describe().coverage['25%']
q1.to_csv("first_q_coverage.csv")

df2=pd.DataFrame()
#for file2 in glob.glob(args.coveragefile+"/*/*"):
for file2 in os.listdir("."):
    if file2.endswith("depth.txt"):
        if os.stat(file2).st_size > 0:
            f=open(file2)
            df2=df2.append(pd.DataFrame(pd.read_table(f,header=None)))
df2.columns=['#rname','position','depth']
df2=df2[df2.depth != 0]
q2=df2.groupby('#rname')[['depth']].describe().depth['25%']
q2.to_csv("first_q_depth.csv")


cutcov=pd.read_csv("first_q_coverage.csv")
#for file1 in glob.glob(args.coveragefile+"/*/*coverage.txt"):
for file1 in os.listdir("."):
    if file1.endswith("coverageall.txt"):
        f=pd.read_csv(file1,sep="\t")
        for i in range(0,len(f)):
            for j in range(0,len(cutcov)):
                if (f['#rname'][i] == cutcov['#rname'][j]) & (f['coverage'][i] >= cutcov['25%'][j]):
                    with open("list.csv",'a') as ofile:
                        ofile.write(file1.split("/")[0]+"\t"+f['#rname'][i]+"\n")
depth=pd.read_csv("first_q_depth.csv")
#for file1 in glob.glob(args.coveragefile+"/*/*depth"):
for file1 in os.listdir("."):
    if file1.endswith("depth.txt"):    
        if os.stat(file1).st_size > 0:
            f=pd.read_csv(file1,sep="\t",header=None)
            f.columns=['#rname','position','depth']
            for i in range(0,len(f)):
                for j in range(0,len(depth)):
            #print(f.head())
                    if (f['#rname'][i] == depth['#rname'][j]) & (f['depth'][i] >= depth['25%'][j]):
                        with open("list1.csv",'a') as ofile1:
                            ofile1.write(file1.split("/")[0]+"\t"+f['#rname'][i]+"\t"+str(f['position'][i])+"\t"+str(f['depth'][i])+"\n")