# vim: set ft=python:

# RNAseq workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to have a STAR index already available within the data/00ref directory
# Also, it is necessary to configure the config.yml file accordingly to include all metadatata required.

import os
import glob

configfile: "config/config.yml"
samples = config["samples"].keys()
genome = config["reference"]["genome_file"]
annotation = config["reference"]["ensembl_gtf"]
starOverhang = config["star_db"]["sjdbOverhang"]

# Build directory structure
paths = ["data/00ref", "data/00reads", "data/00adapters",
        "results/01trim", "results/02abundant","results/03map_reads",
        "results/04dedup","results/05correlation", "results/05bigwig",
        "results/05counts","results/06fastqc_raw",
        "results/06fastqc_trim","results/07multiqc"
        ]
for path in paths:
    os.makedirs(path,exist_ok = True)

# Functions
def get_fq1(wildcards):
            return [ my_files for my_files in glob.glob(f"data/00reads/{wildcards.sample}_*_R1*.fastq.gz")]

def get_fq2(wildcards):
            return [ my_files for my_files in glob.glob(f"data/00reads/{wildcards.sample}_*_R2*.fastq.gz")]

# Set what rules to run locally
localrules: all,
            build_abundant_db,
            fastqc,
            multiqc

rule all:
    # IMPORTANT: output file fo all rule has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input:  #o10 = expand("results/01trim/{s}.{e}.fastq.gz", e=["1P","2P","1U","2U"], s=samples)
            o0 = "results/05counts/read_counts",
            o1 = expand("results/06fastqc_raw/{s}.R1_fastqc.html", s=samples),
            o2 = expand("results/06fastqc_raw/{s}.R2_fastqc.html", s=samples),
            o3 = expand("results/06fastqc_trim/{s}.1P_fastqc.html", s=samples),
            o4 = expand("results/06fastqc_trim/{s}.2P_fastqc.html", s=samples),
            o5 = expand("results/06fastqc_trim/{s}.1U_fastqc.html", s=samples),
            o6 = expand("results/06fastqc_trim/{s}.2U_fastqc.html", s=samples),
            o7 = "results/07multiqc/multiqc_done.flag",
            o8 = "results/05correlation/multiBamSummary.results.npz",
            o9 = expand("results/05bigwig/{s}.bw", s=samples),
            o11 = "data/00ref/SA"
            #expand("00map_reads/{s}.", s=samples) #,
            #expand("00abundant/{s}.fastq.1.gz", s=samples),
            #expand("00abundant/{s}.fastq.2.gz", s=samples)

rule merge_rep:
    input: fq1 = get_fq1,
           fq2 = get_fq2
    output: merge1 = temp("results/00merged_reads/{sample}.R1.fastq.gz"),
            merge2 = temp("results/00merged_reads/{sample}.R2.fastq.gz")
    resources: 
        cpu_per_task = 1,
        partition = "norm",
        time = "3:00:00"
    threads: 1
    shell:
        """
        cat {input.fq1}  > {output.merge1}
        cat {input.fq2}  > {output.merge2}
        """

rule trim:
    input:  fq1 = "results/00merged_reads/{sample}.R1.fastq.gz",
            fq2 = "results/00merged_reads/{sample}.R2.fastq.gz", 

    output:
            fq1P = "results/01trim/{sample}.1P.fastq.gz",
            fq1U = "results/01trim/{sample}.1U.fastq.gz",
            fq2P = "results/01trim/{sample}.2P.fastq.gz",
            fq2U = "results/01trim/{sample}.2U.fastq.gz"
    params: 
            "ILLUMINACLIP:data/00adapters/TruSeq3-PE.fa:2:30:10 MAXINFO:20:0.5 MINLEN:20"
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "3:00:00"
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
    input: config["reference"]["abundant_rna_file"]
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
    
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "3:00:00"
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

rule make_star_index:
    input: gen = f"{genome}", 
        ann = f"{annotation}", 
        so = f"{starOverhang}"
    output: db = "data/00ref/SA",
            db_dir = "data/00ref"
    resources: 
        mem_mb = 32000,
        cpus_per_task = 16,
        partition = "norm",
        time = "3:00:00"
    threads: 16
    shell:
        """
            STAR --runMode genomeGenerate \
             --genomeDir {output.db_dir} \
             --genomeFastaFiles {input.gen} \
             --sjdbGTFfile {input.ann} \
             --sjdbOverhang {input.so}
        """

rule map_reads:
    input: fq1 = "results/02abundant/{sample}.fastq.1.gz",
           fq2 = "results/02abundant/{sample}.fastq.2.gz",
           genome_dir = "data/00ref"
    output: prefix = "results/03map_reads/{sample}",
            bam = "results/03map_reads/{sample}Aligned.sortedByCoord.out.bam"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "3:00:00",
        mem_mb = 32000
    benchmark:
        "benchmarks/map_reads/{sample}.tsv"
    shell:
        """
        STAR --runMode alignReads \
                --runThreadN {threads} \
                --genomeDir {input.genome_dir} \
                --alignSJDBoverhangMin 1 \
                --alignSJoverhangMin 5 \
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
    resources:
        cpus_per_task = 1,
        mem_mb = 4000,
        partition = "norm",
        time = "3:00:00"
    shell:
        """
        picard -Xmx4g MarkDuplicates \
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
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "3:00:00"
    shell:
        """
        bamCoverage -p {threads} -b {input} -o {output} {params}
        """

rule run_correlation_analysis:
    input: expand("results/04dedup/{s}.sorted.dedup.bam", s=samples)
    output: "results/05correlation/multiBamSummary.results.npz"
    benchmark:
        "benchmarks/run_correlation_analysis/correlation.tsv"
    threads: 1
    resources:
        cpus_per_task = 1,
        partition = "norm",
        time = "3:00:00"
    shell: 
        """
        multiBamSummary bins --bamfiles {input} -o {output}
        plotCorrelation -in {output} --whatToPlot heatmap \
         --corMethod pearson -o results/05correlation/heatmap.png \
         --outFileCorMatrix outFileCorMatrix.txt
        """

rule counts:
    input: 
        genome = config["reference"]["genome_file"],
        gtf = config["reference"]["ensembl_gtf"],
        bam = expand("results/04dedup/{s}.sorted.dedup.bam", s=samples)
    output: counts = "results/05counts/read_counts",
            summary = "results/05counts/read_counts.summary"
    params: "-t CDS -g gene_id -O -s 1 -J -R BAM -p --ignoreDup -M --fraction"  # Current params ignore multimappers and duplicated reads
                                                                   # -p = count fragments instead of individual reads
                                                                   # -M = include multi-mapping reads -O count reads mapping overlapping features
                                                                   # --fraction = multimapped reads will be caused as a fraction 
                                                                   #              instead of 1 (1/x where x = numb alignments reported for same read)
    benchmark:
        "benchmarks/counts/counts.tsv"
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "3:00:00"
    shell:
        """
        featureCounts {params} -G {input.genome} -T {threads}\
         -a {input.gtf} \
         -o {output.counts} {input.bam}
        """

rule fastqc:
    input: raw1 = "results/00merged_reads/{sample}.R1.fastq.gz",
           raw2 = "results/00merged_reads/{sample}.R2.fastq.gz",
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

