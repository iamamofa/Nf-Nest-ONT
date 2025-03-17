#!/usr/bin/python3
import os
import sys
import glob
import time
import shutil
import logging
import argparse
import subprocess
import pandas as pd
from pathlib import Path
from itertools import repeat
from multiprocessing import Pool
from nest.bbduk import QualCheck
from nest.alignment import Bwa
from nest.alignment import Bowtie
from nest.alignment import BBMap
from nest.alignment import Snap
from nest.samtools import Samtools
from nest.gatk import GenAnTK
from nest.gatk import Picard
from nest.gatk import FreeBayes
from nest.kestrel import KestrelVar
#from nest.annotater import Annotate
from nest.kestrel import kes_runner
from nest.summarize import Summary
from nest.prepinputs import Prepper
from nest.parsers.vcfReader import Reader 
from nest.parsers.vcfmerge import Merge
from nest.parsers.vcfannotate import Annotate
from nest.parsers.vcfwriter import Writer
import os

def main(r,b,o,v1,v2,v3,m,voi,n):
    ref_path=r
    bed_path=b
    out_path=o
    vcf_path1=v1
    vcf_path2=v2
    vcf_path3=v3
    bam_path=m
    voi_path=voi
    sam_name=n
    main_logger = logging.getLogger('NeST.{0}'.format(sam_name))
    #Filer  and annotate variant calls
    main_logger.debug('Annotating variants')
    ###if file empty don't include in dictionary....}
    #print(vcf_path1)
    #print(os.stat(vcf_path1).st_size)
    #print(os.stat(vcf_path2).st_size)
    #print(os.stat(vcf_path3).st_size)
    vcf_dict={}
    if os.stat(vcf_path1).st_size != 0:
        vcf_dict[vcf_path1]='Samtools'
    if os.stat(vcf_path2).st_size != 0:
        vcf_dict[vcf_path2]='GATK'
    if os.stat(vcf_path3).st_size != 0:
        vcf_dict[vcf_path3]='Freebayes'
    #print(vcf_dict)    
    #vcf_dict = {vcf_path2: 'GATK', vcf_path1: 'Samtools', vcf_path3: 'Freebayes'}
    merger = Merge(out_path, vcf_dict, ref_path)
    merged_vcf = merger.splitter(list(vcf_dict.keys()))[0]
    final_vcf= 'final_{1}.vcf'.format(out_path, sam_name)
    os.rename(merged_vcf, final_vcf)
    return(final_vcf, 0)

if __name__ == '__main__':
    #Define deffault paths and aligner informations
    inp_path = 'Macintosh HD⁩/Users⁩/adminuser⁩/Desktop⁩/CDC⁩/nextflow-NeST2⁩/fq/MaRS_test2'
    ref_def_path = "Macintosh HD⁩/Users⁩/adminuser⁩/Desktop⁩/CDC⁩/nextflow-NeST2⁩/ref⁩/pfalciparum⁩/mdr.fa"
    bbduk_def = 'Macintosh HD/Users⁩/adminuser⁩/anaconda3⁩/bin⁩/bbduk.sh' #"{0}/bbmap/bbduk.sh".format(def_path)
    bbmap_def = 'Macintosh HD/Users⁩/adminuser⁩/anaconda3⁩/bin⁩/bbmap.sh' #"{0}/bbmap/bbmap.sh".format(def_path)
    bwa_def = 'Macintosh HD/Users⁩/adminuser⁩/anaconda3⁩/bin⁩/bwa' #"{0}/bwa/bwa".format(def_path)
    bowtie_def = 'Macintosh HD/Users⁩/adminuser⁩/anaconda3⁩/bin⁩/bowtie2' #"{0}/bowtie2/bowtie2".format(def_path)
    snap_def = 'Macintosh HD⁩/Users⁩/adminuser⁩/anaconda3⁩/pkgs⁩/snap-2013_11_29-0⁩/bin⁩/snap' #"{0}/snap/snap-aligner".format(def_path)
    smt_def = 'Macintosh HD/Users⁩/adminuser⁩/anaconda3⁩/bin⁩/samtools' #"{0}/samtools/samtools".format(def_path)
    bft_def = 'Macintosh HD/Users⁩/adminuser⁩/anaconda3⁩/bin⁩/bcftools' #"{0}/bcftools/bcftools".format(def_path)
    gatk_def = 'Macintosh HD⁩/Users⁩/adminuser⁩/Desktop⁩/CDC⁩/nextflow-NeST2⁩/GenomeAnalysisTK-1.0.5974/GenomeAnalysisTK.jar' #"{0}/GenomeAnalysisTK.jar".format(def_path)
    pic_def = 'Macintosh HD⁩/Users⁩/adminuser⁩/Desktop⁩/CDC⁩/nextflow-NeST2⁩/picard.jar' #"{0}/picard.jar".format(def_path)
    sra_def = 'Macintosh HD⁩/Users⁩/adminuser⁩/Desktop⁩/CDC⁩/Bayesian⁩/lib⁩/sratoolkit⁩/bin⁩/fastq-dump' #'{0}/sratoolkit/bin/fastq-dump'.format(def_path)
    voi_def = "Macintosh HD⁩/Users⁩/adminuser⁩/Desktop⁩/CDC⁩/nextflow-NeST2⁩/ref/pfalciparum/Reportable_SNPs.csv" #'{0}/Reportable_SNPs.csv'.format(ref_def_path)
    #if 'java version "1.8.' in str(subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT).decode('UTF-8').split('\n')[0]):
    java_def = 'java'
    adp_path = "Macintosh HD⁩/Users⁩/adminuser⁩/Desktop⁩/CDC⁩/nextflow-NeST2⁩/ref⁩/pfalciparum⁩/adapters.fa"
    #else:
    #    java_def = "{0}/jdk/bin/java".format(def_path)
    aligner_def = {'bwa' : bwa_def, 'snap' : snap_def, 'bowtie2': bowtie_def, 'bbmap': bbmap_def}
    #Get arguments
    parser = argparse.ArgumentParser(prog='NeST')
    parser.add_argument('-r', '--ref', dest='ref_path', type=str,
                        help='Path to Reference fasta file', required=True)
    parser.add_argument('-b', '--bed', dest='bed_path', type=str,
                        help='Path to Bed file for MDR regions', required=True)
    parser.add_argument('-o', '--outpath', dest='out_path', type=str,
                        help='Path where all outputs will be stored', required=True)
    parser.add_argument('-v1', '--vcf1', dest='vcf_path1', type=str,
                        help='Sample name', required=True)
    parser.add_argument('-v2', '--vcf2', dest='vcf_path2', type=str,
                        help='Sample name', required=True)
    parser.add_argument('-v3', '--vcf3', dest='vcf_path3', type=str,
                        help='Sample name', required=True)
    parser.add_argument('-m', '--bam', dest='bam_path', type=str, required=True,
                         help='The aligner to used by MARs')    
    parser.add_argument('-voi', '--voi', dest='voi_path', type=str, required=True,
                         help='The aligner to used by MARs')    
    parser.add_argument('-name', '--sam_name', dest='sam_name', type=str,
                        help='Sample name', required=True)
    parser.add_argument('--threads', dest='threads', type=int, default=5,
                        help='Number of threads')
    parser.add_argument('--verbose', action='store_true', 
                        help='Increase verbosity of log file')                        
    parser.add_argument('--purge', action='store_true', 
                        help='Remove intermiediate Fastq and alignment files') 

    args = parser.parse_args()
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)
    main(args.ref_path,args.bed_path,args.out_path,args.vcf_path1,args.vcf_path2,args.vcf_path3,args.bam_path,args.voi_path,args.sam_name)