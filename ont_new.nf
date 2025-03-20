#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Path definitions
params.genome = "$baseDir/ref/pfalciparum/New_6_genes.fa"
params.bed = "$baseDir/ref/pfalciparum/New_6_genes.bed"
params.porechop = "porechop"
params.minimap2 = "minimap2"
params.samtools = "samtools"
params.bcftools = "bcftools"
params.medaka = "medaka_consensus"
mode = "$params.input.mode"
vcfmode = "$params.input.vcfmode"
logpath = "$baseDir/.nextflow.log"
workpath = "$baseDir/work"
CTexcelpath = "$baseDir/ANG_2019_TES_master.xlsx"
datasettype = "$params.input.datasettype"
configfiles = "$baseDir/new6genes"
pyscripts = "$baseDir/pyscripts"

// Channel Definitions
fastq_path = Channel.fromFilePairs("${params.input.fastq_path}/barcode*/", size: -1)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
genome_channel = Channel.fromPath(params.genome)


// Workflow Block
workflow {
    // Define input channels
    barcode_folders = Channel.fromPath("${params.input.fastq_path}/*", type: 'dir')

    // Process execution
    concatenateFastq(barcode_folders)
    trimAdapters(concatenateFastq.out)
    qualityControl(trimAdapters.out)
    alignReads(qualityControl.out, genome_channel)
    callVariants(alignReads.out, genome_channel)
    processAlignments(callVariants.out, genome_channel)
    getcoverage(processAlignments.out)
    filterLowCoverage(getcoverage.out)
    generateFinalReport(
        filterLowCoverage.out[0],  // coverage_file
        filterLowCoverage.out[1],  // depth_file
        filterLowCoverage.out[2]   // filtered_regions
    )
    generateVCF(processAlignments.out, genome_channel)
    annotate(generateVCF.out, configfiles)
    vartypes(annotate.out)
    merge(vartypes.out, processAlignments.out, pyscripts)
    filter(merge.out)
    extract(filter.out)
    spread(extract.out, pyscripts)
    snpfilter(spread.out, merge.out, processAlignments.out, pyscripts)
    tojson(snpfilter.out, pyscripts)

    // Pass the collected coverage and depth channels correctly
    finalcoverage = getcoverage.out.collect { it[0] }
    depth = getcoverage.out.collect { it[1] }

    // Pass finalcoverage and depth to getcutoff
    getcutoff(finalcoverage, depth, pyscripts)
    viz(getcutoff.out)
    errortrack(viz.out, workpath, CTexcelpath)
}

// Process Definitions

process concatenateFastq {
    tag { barcode }
    publishDir "${params.output.folder}/${barcode}/concatenated", mode: "copy"

    input:
        path barcode_folder

    output:
        tuple val(barcode), path("${barcode}.fastq")

    script:
        barcode = barcode_folder.name
        """
        cat ${barcode_folder}/*.fastq > ${barcode}.fastq
        """
}

process trimAdapters {
    tag { barcode }
    publishDir "${params.output.folder}/trimmed/${barcode}", pattern: "*.fastq", mode: "copy"

    input:
        tuple val(barcode), path(fastq_group)

    output:
        tuple val(barcode), path("${barcode}_trimmed.fastq")

    script:
        """
        ${params.porechop} -i ${fastq_group} -o ${barcode}_trimmed.fastq --discard_middle
        """
}

process qualityControl {
    tag { barcode }
    publishDir "${params.output.folder}/${barcode}/qc", mode: "copy"

    input:
        tuple val(barcode), path(fastq)

    output:
        tuple val(barcode), path("${barcode}_qc.html")

    script:
        """
        NanoPlot --fastq ${fastq} -o . --html
        mv NanoPlot-report.html ${barcode}_qc.html
        """
}

process alignReads {
    tag { barcode }
    publishDir "${params.output.folder}/aligned/${barcode}", pattern: "*.bam", mode: "copy"
    publishDir "${params.output.folder}/alignedReads/${barcode}/Stats", pattern: "*.txt", mode: "copy"

    input:
        tuple val(barcode), path(fastq)
        path genome

    output:
        tuple val(barcode), path("${barcode}.bam")

    script:
        """
        ${params.minimap2} -ax map-ont ${genome} ${fastq} | ${params.samtools} view -bS - | ${params.samtools} sort -o ${barcode}.bam
        ${params.samtools} index ${barcode}.bam
        """
}

process callVariants {
    tag { barcode }
    publishDir "${params.output.folder}/variants/${barcode}", pattern: "*.vcf", mode: "copy"

    input:
        tuple val(barcode), path(bam)
        path genome

    output:
        tuple val(barcode), path("${barcode}.vcf")

    script:
        """
        ${params.medaka} -i ${bam} -r ${genome} -o ${barcode}_medaka
        ${params.bcftools} mpileup -f ${genome} ${bam} | ${params.bcftools} call -mv -Ov -o ${barcode}.vcf
        """
}

process processAlignments {
    tag { barcode }
    publishDir "${params.output.folder}/alignedReads/${barcode}", pattern: "*.ba*", mode: "copy"
    errorStrategy 'retry'

    input:
        tuple val(barcode), path(bam)
        path genome

    output:
        tuple val(barcode), path("${barcode}_SR.bam")

    script:
        """
        # Add read groups using samtools (ONT-specific)
        ${params.samtools} addreplacerg -r "@RG\\tID:${barcode}\\tSM:${barcode}\\tPL:ONT\\tLB:ONT" -o ${barcode}_withRG.bam ${bam}

        # Sort and index the BAM file
        ${params.samtools} sort -o ${barcode}_SR.bam ${barcode}_withRG.bam
        ${params.samtools} index ${barcode}_SR.bam
        """
}

process getcoverage {
    tag { barcode }
    publishDir "${params.output.folder}/samtoolcoverage/${barcode}", mode: "copy"

    input:
        tuple val(barcode), path(bam)

    output:
        file "${barcode}_coverage.txt"
        file "${barcode}_depth.txt"

    script:
        """
        samtools index ${bam}
        samtools coverage ${bam} -o ${barcode}_coverage.txt
        samtools depth ${bam} -a -o ${barcode}_depth.txt
        """
}

process filterLowCoverage {
    tag { barcode }
    publishDir "${params.output.folder}/filtered_regions/${barcode}", mode: "copy"

    input:
        file coverage_file
        file depth_file

    output:
        file "${barcode}_filtered_regions.bed"
        file "${barcode}_coverage.txt"  // Add this line
        file "${barcode}_depth.txt"     // Add this line

    script:
        """
        # Filter regions with coverage below a threshold (e.g., 10x)
        awk '\$7 < 10' ${coverage_file} > ${barcode}_low_coverage_regions.txt
        awk '\$3 < 10' ${depth_file} > ${barcode}_low_depth_regions.txt
        cat ${barcode}_low_coverage_regions.txt ${barcode}_low_depth_regions.txt | sort -k1,1 -k2,2n | bedtools merge -i - > ${barcode}_filtered_regions.bed
        cp ${coverage_file} ${barcode}_coverage.txt  // Add this line
        cp ${depth_file} ${barcode}_depth.txt        // Add this line
        """
}

process generateFinalReport {
    tag { barcode }
    publishDir "${params.output.folder}/final_reports/${barcode}", mode: "copy"

    input:
        file coverage_file
        file depth_file
        file filtered_regions

    output:
        file "${barcode}_final_report.txt"

    script:
        """
        echo "Coverage Statistics:" > ${barcode}_final_report.txt
        cat ${coverage_file} >> ${barcode}_final_report.txt
        echo "\nDepth Statistics:" >> ${barcode}_final_report.txt
        cat ${depth_file} >> ${barcode}_final_report.txt
        echo "\nFiltered Regions (Low Coverage):" >> ${barcode}_final_report.txt
        cat ${filtered_regions} >> ${barcode}_final_report.txt
        """
}

process generateVCF {
    tag { barcode }
    publishDir "${params.output.folder}/vcf_files/${barcode}", mode: "copy"

    input:
        tuple val(barcode), path(bam)
        path genome

    output:
        file "${barcode}_variants.vcf"

    script:
        """
        ${params.bcftools} mpileup -f ${genome} ${bam} | ${params.bcftools} call -mv -Ov -o ${barcode}_variants.vcf
        """
}

process annotate {
    tag { barcode }
    publishDir "${params.output.folder}/annotate/${barcode}", mode: "copy"
    errorStrategy 'ignore'

    input:
        tuple val(barcode), path(vcf)
        path configpath

    output:
        tuple val(barcode), path("${barcode}.ann.vcf")

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config New_6_genes -interval $baseDir/ref/pfalciparum/annotation.bed ${vcf} > ${barcode}.ann.vcf
        """
}

process vartypes {
    tag { barcode }
    publishDir "${params.output.folder}/vartypes/${barcode}", mode: "copy"

    input:
        tuple val(barcode), path(vcf)

    output:
        tuple val(barcode), path("${barcode}_vartypes.txt")

    script:
        """
        python3 ${pyscripts}/annotate_vartypes.py ${vcf} > ${barcode}_vartypes.txt
        """
}

process merge {
    tag { barcode }
    publishDir "${params.output.folder}/final_vcf/${barcode}", mode: "copy"
    errorStrategy 'ignore'

    input:
        tuple val(barcode), path(vartype_vcf)
        tuple val(barcode), path(bam)
        path pyscripts_path

    output:
        tuple val(barcode), path("final_${barcode}.vcf")

    script:
        """
        samtools index ${bam}
        python3 ${pyscripts_path}/annotate2.py -r $baseDir/ref/pfalciparum/New_6_genes.fa -b $baseDir/ref/pfalciparum/New_6_genes.bed -o ${barcode} -v1 ${vartype_vcf} -m ${bam} -voi $baseDir/ref/pfalciparum/voinew2.csv -name ${barcode}
        """
}

process filter {
    tag { barcode }
    publishDir "${params.output.folder}/filter/${barcode}", mode: "copy"

    input:
        tuple val(barcode), path(final_vcf)

    output:
        tuple val(barcode), path("${barcode}_filtered.vcf")

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar filter -f ${final_vcf} " ( VARTYPE = 'SNP' ) " > ${barcode}_filtered.vcf
        """
}

process extract {
    tag { barcode }
    publishDir "${params.output.folder}/extract/${barcode}", mode: "copy"
    errorStrategy 'ignore'

    input:
        tuple val(barcode), path(filtered_vcf)

    output:
        tuple val(barcode), path("final_${barcode}_ext.vcf")

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar extractFields ${filtered_vcf} CHROM POS REF ALT VARTYPE Confidence Sources DP4 DP AD GEN[*].AD "ANN[*].EFFECT" "ANN[*].HGVS_C" "ANN[*].HGVS_P" > final_${barcode}_ext.vcf
        """
}

process spread {
    tag { barcode }
    publishDir "${params.output.folder}/spread/${barcode}", mode: "copy"

    input:
        tuple val(barcode), path(extracted_vcf)
        path pyscripts_path

    output:
        tuple val(barcode), path("fixedPOS_${barcode}.vcf")

    script:
        """
        python3 ${pyscripts_path}/vcfcsv3.py -n ${extracted_vcf}
        mv ${extracted_vcf}_fixedPOS.vcf fixedPOS_${barcode}.vcf
        """
}

process snpfilter {
    tag { barcode }
    publishDir "${params.output.folder}/snpfilter/${barcode}", pattern: "*.csv", mode: "copy"
    publishDir "${params.output.folder}/log/", pattern: "*.txt", mode: "copy"
    publishDir "${params.output.folder}/forjson/", pattern: "*.csv", mode: "copy"

    input:
        tuple val(barcode), path(spread_vcf)
        tuple val(barcode), path(merged_vcf)
        tuple val(barcode), path(bam)
        path pyscripts_path

    output:
        tuple val(barcode), path("${barcode}.csv")

    script:
        """
        python3 ${pyscripts_path}/snpreport1-54-VAF_SP.py -v1 ${merged_vcf} -v2 ${spread_vcf} -b1 ${bam} -b2 $baseDir/ref/pfalciparum/mdr.bed -o1 ${barcode} -e1 $baseDir/ref/pfalciparum/candidates.xlsx -e2 $baseDir/ref/pfalciparum/voinew4.csv -f1 $baseDir/ref/pfalciparum/New_6_genes.fa
        """
}

process tojson {
    tag { barcode }
    publishDir "${params.output.folder}/combinedjson/", pattern: "*.json", mode: "copy"

    input:
        tuple val(barcode), path(csv_file)
        path pyscripts_path

    output:
        tuple val(barcode), path("wholecombine.json")

    script:
        """
        python3 ${pyscripts_path}/tojson2.py -d1 $baseDir/${params.output.folder}/forjson
        """
}

process getcutoff {
    tag { barcode }
    publishDir "${params.output.folder}/cutoff/", mode: "copy"
    errorStrategy 'ignore'

    input:
        file coverage_files from finalcoverage.collect()
        file depth_files from depth.collect()
        path pyscripts_path

    output:
        file("first_q_coverage.csv")
        file("first_q_depth.csv")
        file("list.csv")
        file("list1.csv")

    script:
        """
        python3 ${pyscripts_path}/getcutoff.py -c $baseDir/${params.output.folder}/samtoolcoverage
        """
}

process viz {
    tag { barcode }
    publishDir "${params.output.folder}/visualization/", mode: "copy"
    errorStrategy 'ignore'

    input:
        path logfile
        path pyscripts_path
        file("*") from snpfilter_out2.collect()
        file list
        file list1

    output:
        file("summary_reportable.csv")
        file("NCBI_feature_column.csv")
        file("PowerBI_input.csv")
        file("cutoff_SNP.csv")
        file("summary.csv")
        file("haplotype.csv")
        file("nextflowlog.txt")

    script:
        """
        python3 ${pyscripts_path}/visualization_final.py -voi $baseDir/ref/pfalciparum/voinew4.csv -cod $baseDir/${params.output.folder}/cutoff/
        mv ${logfile} "nextflowlog.txt"
        """
}

process errortrack {
    tag { barcode }
    publishDir "${params.output.folder}/errortrack/", mode: "copy"
    errorStrategy 'ignore'

    input:
        file logfile
        path work
        path pyscripts_path
        path CT

    output:
        file("out.csv")

    script:
        if (datasettype == 'Angola') {
            """
            python3 ${pyscripts_path}/errortrackn6.py -n1 ${logfile} -w1 ${work} -x1 ${CT}
            """
        } else {
            """
            python3 ${pyscripts_path}/errortrackn7.py -n1 ${logfile} -w1 ${work}
            """
        }
}