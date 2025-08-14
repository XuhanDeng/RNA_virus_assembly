#!/usr/bin/env python3
"""
Yangtze River RNA Virus Pipeline - Snakemake Workflow
RNA virus bioinformatics pipeline for processing RNA NGS sequencing data
"""

# Configuration
configfile: "config/config.yaml"


# Final output rule
rule all:
    input:
        # Step 1: Quality control
        expand("results/1_fastp/{sample}/{sample}_1P.fq.gz", sample=config["samples"]),
        expand("results/1_fastp/{sample}/{sample}_2P.fq.gz", sample=config["samples"]),
        
        # Step 2: rRNA removal
        expand("results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz", sample=config["samples"]),
        expand("results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz", sample=config["samples"]),
        
        # Step 3: Assembly
        expand("results/3_spades/{sample}/scaffolds.fasta", sample=config["samples"]),
        expand("results/3_spades/{sample}/contigs.fasta", sample=config["samples"]),
        
        # Step 4: Reformat assemblies
        expand("results/4_reformat/{sample}_reformat.fasta", sample=config["samples"]),
        
        # Step 5: ORF finding and merging
        "results/5_orfinder/merged_orfs.fasta",
        
        # Step 6: Diamond BLASTP
        "results/6_diamond/blastp_results.tsv"





# Rule 1: Quality control and adapter removal with fastp
rule fastp_qc:
    input:
        r1 = "input/{sample}.R1.fq.gz",
        r2 = "input/{sample}.R2.fq.gz"
    output:
        r1_paired="results/1_fastp/{sample}/{sample}_1P.fq.gz",
        r2_paired="results/1_fastp/{sample}/{sample}_2P.fq.gz",
        r1_unpaired="results/1_fastp/{sample}/{sample}_U1.fq.gz",
        r2_unpaired="results/1_fastp/{sample}/{sample}_U2.fq.gz",
        html="results/1_fastp/{sample}/{sample}.fastp.html",
        json="results/1_fastp/{sample}/{sample}.fastp.json"     
    conda:
        "envs/fastp.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["fastp"]["runtime"],
        cpus_per_task = config["fastp"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/1_fastp/{sample}.log",
        err="log/1_fastp/{sample}.err"
    shell:
        """
        mkdir -p results/1_fastp/{wildcards.sample}
        fastp --thread {threads} \
              --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1_paired} --out2 {output.r2_paired} \
              --unpaired1 {output.r1_unpaired} --unpaired2 {output.r2_unpaired} \
              -h {output.html} -j {output.json} \
              --trim_poly_g --trim_poly_x \
              --qualified_quality_phred {config[fastp][quality_threshold]} \
              --length_required {config[fastp][length_required]} \
              --dont_overwrite > {log.out} 1> {log.err}  # Redirect stderr to log file
        """

# Rule 2: Remove rRNA sequences with ribodetector
rule ribodetector_rrna_removal:
    input:
        r1 = "results/1_fastp/{sample}/{sample}_1P.fq.gz",
        r2 = "results/1_fastp/{sample}/{sample}_2P.fq.gz"
    output:
        r1_nonrrna = "results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz",
        r2_nonrrna = "results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz"
    conda:
        "envs/ribodetector.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["ribodetector"]["runtime"],
        cpus_per_task = config["ribodetector"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/2_ribodetector/{sample}.log",
        err="log/2_ribodetector/{sample}.err"
    shell:
        """
        mkdir -p results/2_ribodetector/{wildcards.sample}
        ribodetector_cpu -t {threads} \
                         -l {config[ribodetector][min_length]} \
                         -i {input.r1} {input.r2} \
                         -e rrna \
                         --chunk_size {config[ribodetector][chunk_size]} \
                         -o {output.r1_nonrrna} {output.r2_nonrrna} > {log.out} 1> {log.err}
        """

# Rule 3: Viral sequence assembly with metaspades
rule spades_assembly:
    input:
        r1 = "results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz",
        r2 = "results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz"
    output:
        scaffolds = "results/3_spades/{sample}/scaffolds.fasta",
        contigs = "results/3_spades/{sample}/contigs.fasta"
    conda:
        "envs/spades.yaml"
    resources:
        mem_mb_per_cpu = config["spades"]["spade_memory"],  # MB
        runtime = config["spades"]["runtime"],
        cpus_per_task = config["spades"]["threads"],
        slurm_partition = config["spades"]["spade_partition"],
        slurm_account = config["account"]
    log:
        out="log/3_spades/{sample}.log",
        err="log/3_spades/{sample}.err"
    shell:
        """
        mkdir -p results/3_spades/{wildcards.sample}
        spades.py --meta \
                  -o results/3_spades/{wildcards.sample} \
                  -1 {input.r1} -2 {input.r2} \
                  -t {threads} -m {resources.mem_mb_per_cpu} \
                  -k {config[spades][kmers]} \
                  --only-assembler > {log.out} 1> {log.err}
        """

# Rule 4: Reformat assembly sequences
rule reformat_assemblies:
    input:
        scaffolds = "results/3_spades/{sample}/scaffolds.fasta"
    output:
        reformatted = "results/4_reformat/{sample}_reformat.fasta"
    conda:
        "envs/seqkit.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = 60,  # minutes
        cpus_per_task = 1,
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/4_reformat/{sample}.log",
        err="log/4_reformat/{sample}.err"
    shell:
        """
        mkdir -p results/4_reformat
        seqkit replace -p "^(.+)$" -r "{wildcards.sample}_scaffold_{{nr}}" \
               --nr-width 10 \
               -o {output.reformatted} \
               {input.scaffolds} > {log.out} 1> {log.err}
        """

#Rule 5 : orfipy to get ORFs for each sample
rule orfipy:
    input:
        fasta = "results/4_reformat/{sample}_reformat.fasta"
    output:
        orfs = "results/5_orfinder/{sample}/{sample}_orfs.fasta"
    conda:
        "envs/orfipy.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],
        runtime = 120,
        cpus_per_task = config["orfipy"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/5_orfinder/{sample}.log",
        err="log/5_orfinder/{sample}.err"
    shell:
        """
        mkdir -p results/5_orfinder/{wildcards.sample}
        orfipy {input.fasta} --dna {output.orfs} --min 10 --max 10000 --procs {threads} --outdir results/5_orfinder/{wildcards.sample}/ > {log.out} 2> {log.err}
        """

#Rule 6: Merge all ORF files
rule merge_orfs:
    input:
        orfs = expand("results/5_orfinder/{sample}/{sample}_orfs.fasta", sample=config["samples"])
    output:
        merged = "results/5_orfinder/merged_orfs.fasta"
    shell:
        """
        cat {input.orfs} > {output.merged}
        """

#Rule 7: Diamond BLASTP against viral RdRP database
rule diamond_blastp:
    input:
        orfs = "results/5_orfinder/merged_orfs.fasta"
    output:
        blastp_results = "results/6_diamond/blastp_results.tsv"
    conda:
        "envs/diamond.yaml"
    params:
        database = config["database"]["diamond_blastp_db"]
    resources:
        mem_mb_per_cpu = config["regular_memory"],
        runtime = config["diamond_blastp"]["runtime"],
        cpus_per_task = config["diamond_blastp"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/6_diamond/blastp.log",
        err="log/6_diamond/blastp.err"
    shell:
        """
        mkdir -p results/6_diamond
        diamond blastp \
            --query {input.orfs} \
            --db {params.database} \
            --out {output.blastp_results} \
            --evalue {config[diamond_blastp][evalue]} \
            --threads {threads} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs 5 > {log.out} 2> {log.err}
        """
