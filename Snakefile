# vim: set ft=python:

# RNAseq workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to have a STAR index already available within the data/00ref directory
# Also, it is necessary to configure the config.yml file accordingly to include all metadatata required.

import os

configfile: "config.yml"
samples = config["samples"].keys()

# Build directory structure
paths = ["data/00ref", "data/00reads", "data/00adapters",
        "results/01trim", "results/02abundant","results/03map_reads",
        "results/04dedup","results/05correlation", "results/05bigwig",
        "results/05counts","results/06fastqc_raw",
        "results/06/fastqc_trim","results/07multiqc"
        ]
for path in paths:
    os.makedirs(path,exist_ok = True)

#localrules: all, clean

rule all:
    # IMPORTANT: output file fo all rule has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input:  
            o0 = "results/05counts/read_counts",
            o1 = expand("results/06fastqc_raw/{s}.R1_fastqc.html", s=samples),
            o2 = expand("results/06fastqc_raw/{s}.R2_fastqc.html", s=samples),
            o3 = expand("results/06fastqc_trim/{s}.1P_fastqc.html", s=samples),
            o4 = expand("results/06fastqc_trim/{s}.2P_fastqc.html", s=samples),
            o5 = expand("results/06fastqc_trim/{s}.1U_fastqc.html", s=samples),
            o6 = expand("results/06fastqc_trim/{s}.2U_fastqc.html", s=samples),
            o7 = "results/07multiqc/multiqc_done.flag",
            o8 = "results/05correlation/multiBamSummary.results.npz",
            o9 = expand("results/05bigwig/{s}.bw", s=samples)
            #expand("00map_reads/{s}.", s=samples) #,
            #expand("00abundant/{s}.fastq.1.gz", s=samples),
            #expand("00abundant/{s}.fastq.2.gz", s=samples)

#def get_fq_files(wildcards):
    

rule merge_rep:
    input: fq1 = "data/00reads/{sample}_L001_R1_001.fastq.gz",
           fq2 = "data/00reads/{sample}_L001_R2_001.fastq.gz"
    output: merge1 = "data/00reads/{sample}.R1.fastq.gz",
            merge2 = "data/00reads/{sample}.R2.fastq.gz"
    shell:
        """
        cat {input.fq1}  > {output.merge1}
        cat {input.fq2}  > {output.merge2}
        """

rule trim:
    input:  fq1 = "data/00reads/{sample}.R1.fastq.gz",
            fq2 = "data/00reads/{sample}.R2.fastq.gz", 

    output:
            fq1P = "results/01trim/{sample}.1P.fastq.gz",
            fq1U = "results/01trim/{sample}.1U.fastq.gz",
            fq2P = "results/01trim/{sample}.2P.fastq.gz",
            fq2U = "results/01trim/{sample}.2U.fastq.gz"
    params: 
            "ILLUMINACLIP:data/00adapters/TruSeq3-PE.fa:2:30:10 MAXINFO:20:0.5 MINLEN:20"
    threads: 16
    log:    "results/01trim/{sample}.log"
    benchmark:
            "benchmarks/trim/{sample}.tsv"
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 \
            {input.fq1} {input.fq2} \
            {output.fq1P} {output.fq1U} {output.fq2P} {output.fq2U} \
            {params} 2> {log}  
        """

rule build_abundant_db:
    input: "data/00ref/abundant_rna.fasta"
    output: "data/00ref/abundant"
    shell:
        """
        bowtie2-build {input} {output}
        touch {output}
        """

rule rm_abundant_rnas:
    input:  fq1P = "results/01trim/{sample}.1P.fastq.gz",
            fq2P = "results/01trim/{sample}.2P.fastq.gz",
            rnas_idx = "data/00ref/abundant"
    output:
        sam = temp("results/02abundant/{sample}.sam"),
        fq1F = "results/02abundant/{sample}.fastq.1.gz",
        fq2F = "results/02abundant/{sample}.fastq.2.gz"
        #fqOut = "results/02abundant/{sample}.fastq.gz" # test.fastq.1.gz test.fastq.2.gz
    threads: 16
    log: 
        metrics = "results/02abundant/{sample}.metrics.txt",
        logs = "results/02abundant/{sample}.log"
    benchmark:
            "benchmarks/rm_abundant_rnas/{sample}.tsv"
    shell:
        """
        bowtie2 --threads {threads} -L 20 -x {input.rnas_idx} \
         --met-file {log.metrics} \
         --un-conc-gz results/02abundant/{wildcards.sample}.fastq.gz \
         -1 {input.fq1P} -2 {input.fq2P} -S {output.sam} 2> {log.logs}
        """
# rule make_star_index:
# shell:
# STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile hg38.ensGene.gtf --sjdbOverhang 74

rule map_reads:
    input: fq1 = "results/02abundant/{sample}.fastq.1.gz",
           fq2 = "results/02abundant/{sample}.fastq.2.gz",
           genome_dir = "data/00ref"
    output: prefix = "results/03map_reads/{sample}",
            bam = "results/03map_reads/{sample}Aligned.sortedByCoord.out.bam"
    threads: 16
    benchmark:
        "benchmarks/map_reads/{sample}.tsv"
    shell:
        """
        STAR --runMode alignReads \
                --runThreadN {threads} \
                --genomeDir {input.genome_dir} \
                --alignSJDBoverhangMin 1 \
                --alignSJoverhangMin 51 \
                --outFilterMismatchNmax 2 \
                --alignEndsType EndToEnd \
                --readFilesIn {input.fq1} {input.fq2} \
                --readFilesCommand zcat --outFileNamePrefix {output.prefix} \
                --quantMode GeneCounts \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes All
        touch {output}                
        """

rule remove_duplicates:
    input: "results/03map_reads/{sample}Aligned.sortedByCoord.out.bam"
    output: "results/04dedup/{sample}.sorted.dedup.bam"
    params: "READ_NAME_REGEX=null REMOVE_DUPLICATES=true"
    log: "results/04dedup/{sample}.sorted.dedup.metrics.txt"
    benchmark:
        "benchmarks/remove_duplicates/{sample}.tsv"

    shell:
        """
        picard MarkDuplicates \
         I={input} \
         O={output} \
         M={log} \
         {params}
        samtools index {output}
        """
rule make_bigwig:
    input: "results/04dedup/{sample}.sorted.dedup.bam"
    output: "results/05bigwig/{sample}.bw"
    params: "--binSize 10 --normalizeUsing BPM" #  + "--filterRNAstrand [forward/reverse]" to plot strand-specific data
    shell:
        """
        bamCoverage -b {input} -o {output} {params}
        """

rule run_correlation_analysis:
    input: expand("results/04dedup/{s}.sorted.dedup.bam", s=config["samples"].keys())
    output: "results/05correlation/multiBamSummary.results.npz"
    benchmark:
        "benchmarks/run_correlation_analysis/correlation.tsv"
    shell: 
        """
        multiBamSummary bins --bamfiles {input} -o {output}
        plotCorrelation -in {output} --whatToPlot heatmap \
         --corMethod pearson -o results/05correlation/heatmap.png \
         --outFileCorMatrix outFileCorMatrix.txt
        """

rule counts:
    input: 
        genome = "data/" + config["reference"]["genome_file"],
        gtf = "data/" + config["reference"]["ensembl_gtf"],
        bam = expand("results/04dedup/{s}.sorted.dedup.bam", s=config["samples"].keys())
    output: counts = "results/05counts/read_counts",
            summary = "results/05counts/read_counts.summary"
    params: "-t CDS -g gene_id -O -s 1 -J -R BAM -p --ignoreDup "  # Current params ignore multimappers and duplicated reads
                                                                   # -p = count fragments instead of individual reads
                                                                   # -M = include multi-mapping reads -O count reads mapping overlapping features
                                                                   # --fraction = multimapped reads will be caused as a fraction 
                                                                   #              instead of 1 (1/x where x = numb alignments reported for same read)
    benchmark:
        "benchmarks/counts/counts.tsv"
    shell:
        """
        featureCounts {params} -G {input.genome} \
         -a {input.gtf} \
         -o {output.counts} {input.bam}
        """

rule fastqc:
    input: raw1 = "data/00reads/{sample}.R1.fastq.gz",
           raw2 = "data/00reads/{sample}.R2.fastq.gz",
           trim1p = "results/01trim/{sample}.1P.fastq.gz",
           trim1u = "results/01trim/{sample}.1U.fastq.gz",
           trim2p = "results/01trim/{sample}.2P.fastq.gz",
           trim2u = "results/01trim/{sample}.2U.fastq.gz"
    output: 
            o1 = "results/06fastqc_raw/{sample}.R1_fastqc.html",
            o2 = "results/06fastqc_raw/{sample}.R2_fastqc.html",
            o3 = "results/06fastqc_trim/{sample}.1P_fastqc.html",
            o4 = "results/06fastqc_trim/{sample}.1U_fastqc.html",
            o5 = "results/06fastqc_trim/{sample}.2P_fastqc.html",
            o6 = "results/06fastqc_trim/{sample}.2U_fastqc.html"
    shell:
        """
        fastqc -o results/06fastqc_raw {input.raw1} {input.raw2}
        fastqc -o results/06fastqc_trim {input.trim1p} {input.trim1u} \
         {input.trim2p} {input.trim2u}
        """
absolute_path = "/home/lorenziha/Downloads/snakemake-class/workflow/"

rule multiqc:
    input: 
           i1 = "results/01trim",
           i2 = "results/02abundant", 
           i3 = "results/03map_reads",
           i4 = "results/04dedup",
           i5 = "results/05counts",
           i6 = "results/05correlation",
           i7 = "results/06fastqc_raw",
           i8 = "results/06fastqc_trim"
    output: "results/07multiqc/multiqc_done.flag"
    shell:
        """
        multiqc -f -d -o results/07multiqc {input.i1} \
                {input.i7} \
                {input.i8} \
                {input.i2} \
                {input.i3} \
                {input.i4} \
                {input.i5} \
                {input.i6}
        touch {output}
        """

