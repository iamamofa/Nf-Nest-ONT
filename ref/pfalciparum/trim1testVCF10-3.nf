#!/usr/bin/env nextflow
BBDUK = file(params.input.bbduk)
params.reads = params.input.fastq_path+'/SRR*_{1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
/*fastq_path.view()*/
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
bbduk_path = "$params.input.bbduk_def"
params.genome = "$baseDir/ref/pfalciparum/6genes2.fa"
params.bed = "$baseDir/ref⁩/pfalciparum⁩/6genes2.bed"
params.fas = "$baseDir/ref/pfalciparum/adapters.fa"
params.gatk= "$baseDir/gatk-4.1.9.0/gatk"
mode = "$params.input.mode"
vcfmode = "$params.input.vcfmode"

process combineFastq {
    publishDir "$params.output.folder/trimFastq/${pair_id}", pattern: "*.fastq", mode : "copy"
    input:
        set val(pair_id), path(fastq_group) from fastq_path

    
    output:
        tuple val(pair_id), path("${pair_id}_R1.fastq"), path("${pair_id}_R2.fastq") into comb_out

    script:
        """
        cat *1.fastq > ${pair_id}_R1.fastq
        cat *2.fastq > ${pair_id}_R2.fastq
        """
}

comb_out.into{preqc_path; trim_path}

process preFastQC {
    publishDir "$params.output.folder/preFastqQC/${sample}", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two) from preqc_path
    
    output:
        tuple val(sample), path("*") into preqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process trimFastq {
    publishDir "$params.output.folder/trimFastq/${sample}", pattern: "*.fastq", mode : "copy"
    publishDir "$params.output.folder/trimFastq/${sample}/Stats", pattern: "*.txt", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two) from trim_path
        path fas from params.fas
    
    output:
        tuple val(sample), path("${sample}_trimmed_R1.fastq"), path("${sample}_trimmed_R2.fastq"), path("${sample}_*.txt") into trim_out

    script:
        """
        $baseDir/bbduk.sh -Xmx1g k=27 hdist=1 edist=1 ktrim=l in=${read_one} in2=${read_two} \\
        out=${sample}_trimmed_R1.fastq ref=${fas} \\
        qtrim=rl minlength=100 \\
        out2=${sample}_trimmed_R2.fastq stats=${sample}_stats.txt 
        """
}

trim_out.into{align_path; postqc_path}

process postFastQC {
    publishDir "$params.output.folder/postFastqQC/${sample}", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from postqc_path
    
    output:
        tuple val(sample), path("*") into postqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process buildIndex {
    tag "$genome.baseName"
    publishDir "$params.output.folder/Bowtie2Index", mode : "copy"

    input:
    path genome from params.genome

    output:
    file 'genome.index*' into index_ch

    script:
    if( mode == 'Bowtie')
    """
    bowtie2-build ${genome} genome.index
    """
}



process alignReads {
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.sam", mode : "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from align_path
        file index from index_ch
        path genome from params.genome
    
    output:
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out

    script:
        index_base = index[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
        if( mode == 'Bowtie')
            """
            bowtie2 --very-sensitive  --dovetail --met-file ${sample}_bmet.txt -p $task.cpus \\
            -x $baseDir/$params.output.folder/Bowtie2Index/$index_base \\
            -1 ${read_one} -2 ${read_two} -p 4 -S ${sample}.sam 
            """
        else if( mode == 'Bwa' )
            """
            bwa mem -t 4 ${genome} ${read_one} ${read_two} > ${sample}.sam 
            """
        else if( mode == 'BBMap' )
            """
            bbmap ref=${genome} in=${read_one} in2=${read_two} out=${sample}.sam 
            """
        else if( mode == 'Snap' )
            """
            snap paired ${genome} ${read_one} ${read_two} -t 4 -o -sam ${sample}.sam $task.cpus \\
            """
}

process processAlignments {
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.ba*", mode : "copy"
    input:
        set val(sample), path(sam_path), path(align_met) from align_out
        path genome from params.genome
    
    output:
        tuple val(sample), path("${sample}_SR.bam") into postal_out
        tuple val(sample), path("${sample}_SR.bam") into postal_out2
        tuple val(sample), path("${sample}_SR.bam") into postal_out3
        tuple val(sample), path("${sample}_SR.bam") into postal_out4

    script:
        """
        java -jar $baseDir/picard.jar AddOrReplaceReadGroups -I ${sam_path} -O ${sample}_SR.bam -SORT_ORDER coordinate --CREATE_INDEX true \\
        -LB ExomeSeq -DS ExomeSeq -PL Illumina -CN AtlantaGenomeCenter -DT 2016-08-24 -PI null -ID ${sample} \\
        -PG ${sample} -PM ${sample} -SM ${sample} -PU HiSeq2500
        """
    
}

process GenerateVCF {
    publishDir "$params.output.folder/GenerateVCFSam/${sample}", pattern: "*.vcf", mode : "copy"
    input:
        set val(sample), path(bam_path) from postal_out
        path genome from params.genome

    output:
        tuple val(sample), path("${sample}-1.vcf"), path("${sample}-2.vcf"), path("${sample}-3.vcf") into vcf_out1

    script:
        """
        bcftools mpileup -f ${genome} ${bam_path} > ${sample}.mpileup
        bcftools call -vm ${sample}.mpileup > ${sample}-1.vcf
        $baseDir/gatk-4.1.9.0/gatk HaplotypeCaller -R $baseDir/ref/pfalciparum/6genes2.fa -I ${bam_path} -O ${sample}-2.vcf
        freebayes -f ${genome} ${bam_path} > ${sample}-3.vcf
        """
}

/*process GenerateVCFHap {
    publishDir "$params.output.folder/GenerateVCFHap/${sample}", pattern: "*.vcf", mode : "copy"
    input:
        set val(sample), path(bam_path) from postal_out2
        path genome from params.genome
        path gatk from params.gatk

    output:
        tuple val(sample), path("${sample}-2.vcf") into vcf_out2

    script:
        """
        $baseDir/gatk-4.1.9.0/gatk HaplotypeCaller -R $baseDir/ref/pfalciparum/6genes2.fa -I ${bam_path} -O ${sample}-2.vcf \\
        --output-mode EMIT_VARIANTS_ONLY
        """
}

process GenerateVCFHFree {
    publishDir "$params.output.folder/GenerateVCFFree/${sample}", pattern: "*.vcf", mode : "copy"
    input:
        set val(sample), path(bam_path) from postal_out3
        path genome from params.genome
        path gatk from params.gatk

    output:
        tuple val(sample), path("${sample}-3.vcf") into vcf_out3

    script:
        """
        freebayes -f ${genome} ${bam_path} > ${sample}-3.vcf
        """
}*/

/* process extract {
    publishDir "$params.output.folder/extract/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from vcf_out1
        set val(sample), path(vcf_path2) from vcf_out2
        set val(sample), path(vcf_path3) from vcf_out3

    output:
        tuple val(sample), path("final_${vcf_path1}ext.vcf") into extract_out1
        tuple val(sample), path("final_${vcf_path2}ext.vcf") into extract_out2
        tuple val(sample), path("final_${vcf_path3}ext.vcf") into extract_out3

    script:
    """
    java -jar $baseDir/snpEFF/SnpSift.jar extractFields ${vcf_path1} CHROM POS REF ALT "ANN[*].EFFECT" > final_${vcf_path1}ext.vcf
    java -jar $baseDir/snpEFF/SnpSift.jar extractFields ${vcf_path2} CHROM POS REF ALT "ANN[*].EFFECT" > final_${vcf_path2}ext.vcf 
    java -jar $baseDir/snpEFF/SnpSift.jar extractFields ${vcf_path3} CHROM POS REF ALT "ANN[*].EFFECT" > final_${vcf_path3}ext.vcf
    """
} */


process annotate {
    publishDir "$params.output.folder/annotate/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1), path(vcf_path2), path(vcf_path3) from vcf_out1

    output:
        tuple val(sample), path("${sample}-1.ann.vcf"), path("${sample}-2.ann.vcf"),path("${sample}-3.ann.vcf") into anno_out1 

    script:
        """
        java -Xmx8g -jar $baseDir/snpEff/snpEff.jar 6genes2 ${vcf_path1} > ${sample}-1.ann.vcf
        java -Xmx8g -jar $baseDir/snpEff/snpEff.jar 6genes2 ${vcf_path2} > ${sample}-2.ann.vcf
        java -Xmx8g -jar $baseDir/snpEff/snpEff.jar 6genes2 ${vcf_path3} > ${sample}-3.ann.vcf
        """

}

/*process vartpype {
    publishDir "$params.output.folder/vartype/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from anno_out1
        set val(sample), path(vcf_path2) from anno_out2
        set val(sample), path(vcf_path3) from anno_out3

    output:
        tuple val(sample), path("${sample}_1.vartype.vcf") into vartype_out1
        tuple val(sample), path("${sample}_2.vartype.vcf") into vartype_out2
        tuple val(sample), path("${sample}_3.vartype.vcf") into vartype_out3

    script:
        """
        java -jar $baseDir/snpEFF/SnpSift.jar varType ${vcf_path1} > ${sample}_1.vartype.vcf
        java -jar $baseDir/snpEFF/SnpSift.jar varType ${vcf_path2} > ${sample}_2.vartype.vcf
        java -jar $baseDir/snpEFF/SnpSift.jar varType ${vcf_path3} > ${sample}_3.vartype.vcf
        """
}*/

/* process extract2 {
    publishDir "$params.output.folder/extract/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from anno_out1
        set val(sample), path(vcf_path2) from anno_out2
        set val(sample), path(vcf_path3) from anno_out3

    output:
        tuple val(sample), path("final_${vcf_path1}ext.vcf") into extract_out1
        tuple val(sample), path("final_${vcf_path2}ext.vcf") into extract_out2
        tuple val(sample), path("final_${vcf_path3}ext.vcf") into extract_out3

    script:
    """
    java -jar $baseDir/snpEFF/SnpSift.jar extractFields ${vcf_path1} CHROM POS REF ALT varType "ANN[*].EFFECT" > final_${vcf_path1}ext.vcf
    java -jar $baseDir/snpEFF/SnpSift.jar extractFields ${vcf_path2} CHROM POS REF ALT varType "ANN[*].EFFECT" > final_${vcf_path2}ext.vcf 
    java -jar $baseDir/snpEFF/SnpSift.jar extractFields ${vcf_path3} CHROM POS REF ALT varType "ANN[*].EFFECT" > final_${vcf_path3}ext.vcf
    """
} */

process merge {
    publishDir "$params.output.folder/final_vcf/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1), path(vcf_path2), path(vcf_path3) from anno_out1
        set val(sample), path(bam_path) from postal_out4

    output:
        tuple val(sample), path("final_${sample}.vcf") into merge_out

    script:
        """
        samtools index ${bam_path}
        python $baseDir/annotate.py -r $baseDir/ref/pfalciparum/6genes2.fa -b $baseDir/6genes2.bed -o ${sample} -v1 ${vcf_path1} -v2 ${vcf_path2} -v3 ${vcf_path3} -m ${bam_path} -voi $baseDir/ref/pfalciparum/Reportable_SNPs.csv -name ${sample}
        """
}

process vartpype {
    publishDir "$params.output.folder/vartype/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from merge_out

    output:
        tuple val(sample), path("${sample}_1.vartype.vcf") into vartype_out1

    script:
        """
        java -jar $baseDir/snpEFF/SnpSift.jar varType ${vcf_path1} > ${sample}_1.vartype.vcf
        """
}

process filter {
    publishDir "$params.output.folder/filter/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from vartype_out1

    output:
        tuple val(sample), path("${sample}_filtered.vcf") into filter_out

    script:
        """
        java -jar $baseDir/snpEFF/SnpSift.jar filter -f ${vcf_path1} " ( VARTYPE = 'SNP' ) " > ${sample}_filtered.vcf
        """
}

process extract {
    publishDir "$params.output.folder/extract/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from filter_out

    output:
        tuple val(sample), path("final_${vcf_path1}ext.vcf") into extract_out

    script:
    """
    java -jar $baseDir/snpEFF/SnpSift.jar extractFields ${vcf_path1} CHROM POS REF ALT VARTYPE Confidence Sources AF QD FS MQRankSum MQ "ANN[*].EFFECT" "ANN[*].HGVS_C" "ANN[*].HGVS_P" > final_${vcf_path1}ext.vcf
    """
}

process spread {
    publishDir "$params.output.folder/spread/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path) from extract_out

    output:
        tuple val(sample), path("fixedPOS${vcf_path}") into spread_out

    script:
    """
    python $baseDir/vcfcsv3.py -n ${vcf_path}
    """
}

/*process coverage {
    publishDir "$params.output.folder/coverage/${sample}", mode : "copy"

    input:
        set val(sample), path(vcf_path1) from filter_out

    output:
        tuple val(sample), path("${sample}.txt") into coverage_out

    script:
        """
        bedtools coverage -a $baseDir/ref/pfalciparum/mdr.bed -b ${vcf_path1} > ${sample}.txt
        """
}*/