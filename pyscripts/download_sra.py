import subprocess
import os
import sys
# samples correspond to Het_1, Het_2, Imm_1, Imm_2
sra_numbers = open(sys.argv[1])

# # this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
# for sra_id in sra_numbers:
#     sra_id=sra_id.rstrip()
#     print ("Currently downloading: " + sra_id)
#     prefetch = "prefetch " + sra_id
#     print ("The command used was: " + prefetch)
#     subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    sra_id=sra_id.rstrip()
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
       