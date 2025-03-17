#!/usr/bin/env nextflow
params.reads = params.input.fastq_path + 'barcode*/'  
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
params.genome = "$baseDir/ref/pfalciparum/New_6_genes.fa"
params.bed = "$baseDir/ref⁩/pfalciparum⁩/New_6_genes.bed"
params.porechop = "porechop"
params.minimap2 = "minimap2"
params.samtools = "samtools"
params.bcftools = "bcftools"
params.medaka = "medaka_consensus"
mode = "$params.input.mode"
vcfmode = "$params.input.vcfmode"
logpath="$baseDir/.nextflow.log"
workpath="$baseDir/work"
CTexcelpath="$baseDir/ANG_2019_TES_master.xlsx"
datasettype = "$params.input.datasettype"
configfiles="$baseDir/new6genes"
pyscripts="$baseDir/pyscripts"

// Step 1: Collect barcode folders
barcode_folders = Channel.fromPath("${params.reads}/*", type: 'dir')

// Step 2: Concatenate FASTQ files per barcode folder
process concatenateFastq {
    tag { barcode }
    publishDir "${params.output}/${barcode}/concatenated", mode: "copy"

    input:
        path barcode_folder from barcode_folders

    output:
        tuple val(barcode), path("${barcode}.fastq") into concatenated_fastq

    script:
        barcode = barcode_folder.name
        """
        cat ${barcode_folder}/*.fastq > ${barcode}.fastq
        """
}

process trimAdapters {
    publishDir "$params.output.folder/trimmed/${barcode}", pattern: "*.fastq", mode: "copy"
    
    input:
        set val(barcode), path(fastq_group) from fastq_path

    output:
        tuple val(barcode), path("${barcode}_trimmed.fastq") into trimmed_out

    script:
        """
        ${params.porechop} -i ${fastq_group} -o ${barcode}_trimmed.fastq --discard_middle
        """
}

trimmed_out.into { align_path }


process qualityControl {
    tag { barcode }
    publishDir "${params.output}/${barcode}/qc", mode: "copy"

    input:
        tuple val(barcode), path(fastq) from trimmed_fastq

    output:
        tuple val(barcode), path("${barcode}_qc.html") into qc_out

    script:
        """
        NanoPlot --fastq ${fastq} -o . --html
        mv NanoPlot-report.html ${barcode}_qc.html
        """
}

process alignReads {
    publishDir "$params.output.folder/aligned/${barcode}", pattern: "*.bam", mode: "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    
    input:
        set val(barcode), path(fastq) from align_path
        path genome from params.genome
    
    output:
        tuple val(barcode), path("${barcode}.bam") into aligned_out

    script:
        """
        ${params.minimap2} -ax map-ont ${genome} ${fastq} | ${params.samtools} view -bS - | ${params.samtools} sort -o ${barcode}.bam
        ${params.samtools} index ${barcode}.bam
        """
}

aligned_out.into { variant_call_path }


process callVariants {
    publishDir "$params.output.folder/variants/${barcode}", pattern: "*.vcf", mode: "copy"
    
    input:
        set val(barcode), path(bam) from variant_call_path
        path genome from params.genome
    
    output:
        tuple val(barcode), path("${barcode}.vcf") into vcf_out

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
        set val(barcode), path(bam) from aligned_out
        path genome from params.genome
    
    output:
        tuple val(barcode), path("${barcode}_SR.bam") into postal_out
        tuple val(barcode), path("${barcode}_SR.bam") into postal_out2
        tuple val(barcode), path("${barcode}_SR.bam") into postal_out3
        tuple val(barcode), path("${barcode}_SR.bam") into postal_out4
        tuple val(barcode), path("${barcode}_SR.bam") into postal_out5
        tuple val(barcode), path("${barcode}_SR.bam") into postal_out6

    script:
        """
        # Add read groups using samtools (ONT-specific)
        ${params.samtools} addreplacerg -r "@RG\\tID:${barcode}\\tSM:${barcode}\\tPL:ONT\\tLB:ONT" -o ${barcode}_withRG.bam ${bam}

        # Sort and index the BAM file
        ${params.samtools} sort -o ${barcode}_SR.bam ${barcode}_withRG.bam
        ${params.samtools} index ${barcode}_SR.bam
        """
}

// Output from processAlignments
aligned_out.into { variant_call_path; postal_out; postal_out2; postal_out3; postal_out4; postal_out5; postal_out6 }

// getcoverage process
process getcoverage {
    tag { barcode }
    publishDir "${params.output.folder}/samtoolcoverage/${barcode}", mode: "copy"
    
    input:
    set val(barcode), path(bam) from postal_out6 

    output:
    file "${barcode}_coverage.txt" into coverage
    file "${barcode}_depth.txt" into depth

    script:
    """
    samtools index ${bam}
    samtools coverage ${bam} -o ${barcode}_coverage.txt
    samtools depth ${bam} -a -o ${barcode}_depth.txt
    """
}

// Output from getcoverage
coverage.into { coverage_plots; low_coverage_filter }
depth.into { depth_stats }

// Process: Generate Coverage Plots
process generateCoveragePlots {
    tag { barcode }
    publishDir "${params.output.folder}/coverage_plots/${barcode}", mode: "copy"

    input:
    file coverage_file from coverage_plots

    output:
    file "${barcode}_coverage_plot.png" into coverage_plots_out

    script:
    """
    # Use a plotting tool like matplotlib (Python) or R to generate coverage plots
    python3 $baseDir/pyscripts/generate_coverage_plot.py -i ${coverage_file} -o ${barcode}_coverage_plot.png
    """
}

// Process: Filter Regions with Low Coverage
process filterLowCoverage {
    tag { barcode }
    publishDir "${params.output.folder}/filtered_regions/${barcode}", mode: "copy"

    input:
    file coverage_file from low_coverage_filter
    file depth_file from depth_stats

    output:
    file "${barcode}_filtered_regions.bed" into filtered_regions_out

    script:
    """
    # Filter regions with coverage below a threshold (e.g., 10x)
    awk '\$7 < 10' ${coverage_file} > ${barcode}_low_coverage_regions.txt
    awk '\$3 < 10' ${depth_file} > ${barcode}_low_depth_regions.txt
    cat ${barcode}_low_coverage_regions.txt ${barcode}_low_depth_regions.txt | sort -k1,1 -k2,2n | bedtools merge -i - > ${barcode}_filtered_regions.bed
    """
}

// Process: Integrate Coverage Statistics into a Final Report
process generateFinalReport {
    tag { barcode }
    publishDir "${params.output.folder}/final_reports/${barcode}", mode: "copy"

    input:
    file coverage_file from coverage
    file depth_file from depth
    file filtered_regions from filtered_regions_out

    output:
    file "${barcode}_final_report.txt" into final_reports_out

    script:
    """
    # Combine coverage, depth, and filtered regions into a final report
    echo "Coverage Statistics:" > ${barcode}_final_report.txt
    cat ${coverage_file} >> ${barcode}_final_report.txt
    echo "\nDepth Statistics:" >> ${barcode}_final_report.txt
    cat ${depth_file} >> ${barcode}_final_report.txt
    echo "\nFiltered Regions (Low Coverage):" >> ${barcode}_final_report.txt
    cat ${filtered_regions} >> ${barcode}_final_report.txt
    """
}

// Process: Generate VCF File
process generateVCF {
    tag { barcode }
    publishDir "${params.output.folder}/vcf_files/${barcode}", mode: "copy"

    input:
    set val(barcode), path(bam) from postal_out6
    path genome from params.genome

    output:
    file "${barcode}_variants.vcf" into vcf_out

    script:
    """
    # Call variants using bcftools
    ${params.bcftools} mpileup -f ${genome} ${bam} | ${params.bcftools} call -mv -Ov -o ${barcode}_variants.vcf
    """
}

// Output from generateVCF
generateVCF.out.vcf_out.into { vcf_out }

// annotate process
process annotate {
    tag { barcode }
    publishDir "${params.output.folder}/annotate/${barcode}", mode: "copy"
    errorStrategy 'ignore'

    input:
        set val(barcode), path(vcf) from vcf_out
        path(configpath) from configfiles

    output:
        tuple val(barcode), path("${barcode}.ann.vcf") into anno_out

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/snpEff.jar -c ${configpath}/snpEff.config New_6_genes -interval $baseDir/ref/pfalciparum/annotation.bed ${vcf} > ${barcode}.ann.vcf
        """
}


// Output from annotate
annotate.out.anno_out.into { anno_out }

// vartype process
process vartype {
    tag { barcode }
    publishDir "${params.output.folder}/vartype/${barcode}", mode: "copy"

    input:
        set val(barcode), path(ann_vcf) from anno_out

    output:
        tuple val(barcode), path("${barcode}.vartype.vcf") into vartype_out

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar varType ${ann_vcf} > ${barcode}.vartype.vcf
        """
}

// Output from vartype
vartype.out.vartype_out.into { vartype_out }

// merge process
process merge {
    tag { barcode }
    publishDir "${params.output.folder}/final_vcf/${barcode}", mode: "copy"
    errorStrategy 'ignore'

    input:
        set val(barcode), path(vartype_vcf) from vartype_out
        set val(barcode), path(bam) from postal_out4
        path(pyscripts_path) from pyscripts

    output:
        tuple val(barcode), path("final_${barcode}.vcf") into merge_out
        tuple val(barcode), path("final_${barcode}.vcf") into merge_out2
        tuple val(barcode), path("final_${barcode}.vcf") into merge_out3

    script:
        """
        samtools index ${bam}
        python3 ${pyscripts_path}/annotate2.py -r $baseDir/ref/pfalciparum/New_6_genes.fa -b $baseDir/ref/pfalciparum/New_6_genes.bed -o ${barcode} -v1 ${vartype_vcf} -m ${bam} -voi $baseDir/ref/pfalciparum/voinew2.csv -name ${barcode}
        """
}

// Output from merge
merge.out.merge_out.into { merge_out }

// filter process
process filter {
    tag { barcode }
    publishDir "${params.output.folder}/filter/${barcode}", mode: "copy"

    input:
        set val(barcode), path(final_vcf) from merge_out

    output:
        tuple val(barcode), path("${barcode}_filtered.vcf") into filter_out

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar filter -f ${final_vcf} " ( VARTYPE = 'SNP' ) " > ${barcode}_filtered.vcf
        """
}

// Previous processes (concatenateFastq, trimAdapters, alignReads, processAlignments, getcoverage, generateVCF, annotate, filter, etc.)

// Output from filter
filter.out.filter_out.into { filter_out }

// extract process
process extract {
    tag { barcode }
    publishDir "${params.output.folder}/extract/${barcode}", mode: "copy"
    errorStrategy 'ignore'

    input:
        set val(barcode), path(filtered_vcf) from filter_out

    output:
        tuple val(barcode), path("final_${barcode}_ext.vcf") into extract_out

    script:
        """
        java -Xmx8g -jar $baseDir/tmp/snpEff/SnpSift.jar extractFields ${filtered_vcf} CHROM POS REF ALT VARTYPE Confidence Sources DP4 DP AD GEN[*].AD "ANN[*].EFFECT" "ANN[*].HGVS_C" "ANN[*].HGVS_P" > final_${barcode}_ext.vcf
        """
}

// Output from extract
extract.out.extract_out.into { extract_out }

// spread process
process spread {
    tag { barcode }
    publishDir "${params.output.folder}/spread/${barcode}", mode: "copy"

    input:
        set val(barcode), path(extracted_vcf) from extract_out
        path(pyscripts_path) from pyscripts

    output:
        tuple val(barcode), path("fixedPOS_${barcode}.vcf") into spread_out

    script:
        """
        python3 ${pyscripts_path}/vcfcsv3.py -n ${extracted_vcf}
        mv ${extracted_vcf}_fixedPOS.vcf fixedPOS_${barcode}.vcf
        """
}

// Output from spread
spread.out.spread_out.into { spread_out }

// snpfilter process
process snpfilter {
    tag { barcode }
    publishDir "${params.output.folder}/snpfilter/${barcode}", pattern: "*.csv", mode: "copy"
    publishDir "${params.output.folder}/log/", pattern: "*.txt", mode: "copy"
    publishDir "${params.output.folder}/forjson/", pattern: "*.csv", mode: "copy"

    input:
        set val(barcode), path(spread_vcf) from spread_out
        set val(barcode), path(merged_vcf) from merge_out2
        set val(barcode), path(bam) from postal_out5
        path(pyscripts_path) from pyscripts

    output:
        tuple val(barcode), path("${barcode}.csv") into snpfilter_out
        tuple val(barcode), path("${barcode}.csv") into snpfilter_out2

    script:
        """
        python3 ${pyscripts_path}/snpreport1-54-VAF_SP.py -v1 ${merged_vcf} -v2 ${spread_vcf} -b1 ${bam} -b2 $baseDir/ref/pfalciparum/mdr.bed -o1 ${barcode} -e1 $baseDir/ref/pfalciparum/candidates.xlsx -e2 $baseDir/ref/pfalciparum/voinew4.csv -f1 $baseDir/ref/pfalciparum/New_6_genes.fa
        """
}

// Previous processes (concatenateFastq, trimAdapters, alignReads, processAlignments, getcoverage, generateVCF, annotate, filter, snpfilter, etc.)

// Output from snpfilter
snpfilter.out.snpfilter_out.into { snpfilter_out }
snpfilter.out.snpfilter_out2.into { snpfilter_out2 }

// tojson process
process tojson {
    tag { barcode }
    publishDir "${params.output.folder}/combinedjson/", pattern: "*.json", mode: "copy"

    input:
        set val(barcode), path(csv_file) from snpfilter_out
        path(pyscripts_path) from pyscripts

    output:
        tuple val(barcode), path("wholecombine.json") into tojson_out

    script:
        """
        python3 ${pyscripts_path}/tojson2.py -d1 $baseDir/${params.output.folder}/forjson
        """
}

// Output from tojson
tojson.out.tojson_out.into { tojson_out }

// getcutoff process
process getcutoff {
    tag { barcode }
    publishDir "${params.output.folder}/cutoff/", mode: "copy"
    errorStrategy 'ignore'

    input:
        file("*") from finalcoverage.collect()
        file("*") from depth.collect()
        path(pyscripts_path) from pyscripts

    output:
        file("first_q_coverage.csv")
        file("first_q_depth.csv")
        file("list.csv") into list
        file("list1.csv") into list1

    script:
        """
        python3 ${pyscripts_path}/getcutoff.py -c $baseDir/${params.output.folder}/samtoolcoverage
        """
}

// Output from getcutoff
getcutoff.out.list.into { list }
getcutoff.out.list1.into { list1 }

// viz process
process viz {
    tag { barcode }
    publishDir "${params.output.folder}/visualization/", mode: "copy"
    errorStrategy 'ignore'

    input:
        path(logfile) from logpath
        path(pyscripts_path) from pyscripts
        file("*") from snpfilter_out2.collect()
        file(list) from list
        file(list1) from list1

    output:
        file("summary_reportable.csv")
        file("NCBI_feature_column.csv")
        file("PowerBI_input.csv")
        file("cutoff_SNP.csv")
        file("summary.csv")
        file("haplotype.csv")
        file("nextflowlog.txt") into logpath1

    script:
        """
        python3 ${pyscripts_path}/visualization_final.py -voi $baseDir/ref/pfalciparum/voinew4.csv -cod $baseDir/${params.output.folder}/cutoff/
        mv ${logfile} "nextflowlog.txt"
        """
}

















process errortrack {
    publishDir "$params.output.folder/errortrack/", mode : "copy"
    errorStrategy 'ignore'
    input:

        file(logfile) from logpath1

        path(work) from workpath
        path(pyscripts_path) from pyscripts

        path(CT) from CTexcelpath
    output:
        file("out.csv")
    script:
        if( datasettype == 'Angola')
            """
            python3 ${pyscripts_path}/errortrackn6.py -n1 ${logfile} -w1 ${work} -x1 ${CT}
            """
        else
            """
           python3 ${pyscripts_path}/errortrackn7.py -n1 ${logfile} -w1 ${work}
            """
}

// process sendemail{
//     """
//     sendMail {
//     to 'nej1@cdc.gov'
//     from 'nej1@cdc.gov'
//     attach "$params.output.folder/visualization/PowerBI_input.csv"
//     attach "$params.output.foler/snpfilter/*/*.csv"
//     subject 'Your nfNeST job is done'

//     '''
//     Hi there,
//     Please find the files attached for nfNeST results
//     Thanks,
//     Subin
//     '''
    
//     }
//     """
// }