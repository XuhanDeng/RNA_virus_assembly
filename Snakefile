#!/usr/bin/env python3
"""
Yangtze River RNA Virus Pipeline - Snakemake Workflow
RNA virus bioinformatics pipeline for processing RNA NGS sequencing data
"""

# Configuration
configfile: "config/config.yaml"


# Final output rule
rule all:  
    '''
    Complete pipeline targets (commented out intermediate steps):
    
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



    # Step 5: EsViritu identification
        expand("results/5_esviritu/{sample}/{sample}.detected_virus.info.tsv", sample=config["samples"]),
    # Final outputs only - Snakemake automatically resolves dependencies
        "results/7_merged/esviritu_VIR_table.xlsx"
            '''    
    input:
    # Step 6: Bowtie2 alignment
        expand("results/6_bowtie2/reduntant/2_bowtie2_alignment_bam/{sample}.bam", sample=config["samples"]),
    # Step 7: CoverM analysis
        expand("results/6_bowtie2/reduntant/3_coverm/count/{sample}/{sample}.count.tsv", sample=config["samples"]),
        expand("results/6_bowtie2/reduntant/3_coverm/coverage/{sample}/{sample}.coverage.tsv", sample=config["samples"]),
        expand("results/6_bowtie2/reduntant/3_coverm/tpm/{sample}/{sample}.tpm.tsv", sample=config["samples"]),


# Rule 0: Download EsViritu database (runs only once)
rule download_esviritu_database:
    output:
        db_flag = "databases/esviritu_DB/v3.1.1/database.ready"
    conda:
        "envs/esviritu.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = 60,  # minutes
        cpus_per_task = 1,
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/0_database/esviritu_db_download.log",
        err="log/0_database/esviritu_db_download.err"
    shell:
        """
        # Create directories
        mkdir -p databases/esviritu_DB
        mkdir -p log/0_database
        cd databases/esviritu_DB
        
        # Download database (~400 MB)
        echo "Downloading EsViritu database (v3.1.1, ~400 MB)..." > ../../{log.out}
        wget https://zenodo.org/records/15723755/files/esviritu_db_v3.1.1.tar.gz 2>> ../../{log.err}
        
        # Extract database
        echo "Extracting database..." >> ../../{log.out}
        tar -xvf esviritu_db_v3.1.1.tar.gz >> ../../{log.out} 2>> ../../{log.err}
        
        # Clean up
        echo "Cleaning up..." >> ../../{log.out}
        rm esviritu_db_v3.1.1.tar.gz
        
        # Create flag file to indicate completion
        touch v3.1.1/database.ready
        echo "Database setup complete!" >> ../../{log.out}
        echo "Database location: $(pwd)/v3.1.1" >> ../../{log.out}
        """

'''
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
'''
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

# Rule 5: Read-based viral identification with EsViritu
rule esviritu_identification:
    input:
        r1 = "results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz",
        r2 = "results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz",
        db_flag = "databases/esviritu_DB/v3.1.1/database.ready"
    output:
        results = "results/5_esviritu_merge/{sample}.done"
    conda:
        "envs/esviritu.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["esviritu"]["runtime"],
        cpus_per_task = config["esviritu"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/5_esviritu/{sample}.log",
        err="log/5_esviritu/{sample}.err"
    shell:
        """
        mkdir -p results/5_esviritu_merge
        EsViritu -r {input.r1} {input.r2} \
                 -s {wildcards.sample} \
                 -t {threads} \
                 -o results/5_esviritu_merge \
                 -p {config[esviritu][mode]} > {log.out} 2>> {log.err}
        
        # Create completion flag
        touch {output.results}
        """



# Rule 7: Merge EsViritu results
rule merge_esviritu_results:
    input:
        results = expand("results/5_esviritu_merge/{sample}.done", sample=config["samples"])
    output:
        merged = "results/7_merged/esviritu_VIR_table.xlsx"
    conda:
        "envs/esviritu.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = 30,  # minutes
        cpus_per_task = 1,
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/7_merged/merge_results.log",
        err="log/7_merged/merge_results.err"
    shell:
        """
        mkdir -p results/7_merged \
        summarize_esv_runs {input.results} \
        touch {output.merged}   > {log.out} 2>> {log.err}
        """

#Rule 8: establish bowtie2 reference database
rule bowtie2_database:
    input:
        db_flag = config["databases"]["ncbi_vir_reduntant"]
    output:
        ncbi_vir_reduntant_index = "results/6_bowtie2/reduntant/reduntant_index.done"
    conda:"envs/bowtie2_samtools.yaml"
    log:
        out="log/6_bowtie2/reduntant/bowtie2_database_index_redudant.log",
        err="log/6_bowtie2/reduntant/bowtie2_database_index_redudant.err"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = 60,  # minutes
        cpus_per_task = 8,
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    shell:
        """
        mkdir -p results/6_bowtie2/reduntant/reduntant_index
        bowtie2-build {input.db_flag} \
                     results/6_bowtie2/reduntant/reduntant_index/ncbi_vir_reduntant_index \
                      --threads {resources.cpus_per_task} > {log.out} 2> {log.err}
        touch {output.ncbi_vir_reduntant_index}
        echo "Bowtie2 index for NCBI redundant database created successfully."
        """ 
#Rule 9 :very_sensitive bowtie2 to reduntant
rule bowtie2_very_sensitive:
    input:
        r1 = "results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz",
        r2 = "results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz",
        db_index = "results/6_bowtie2/reduntant/reduntant_index.done"
    output:
        sam = temp("results/6_bowtie2/reduntant/1_bam_file/{sample}.sam"),
        bam = "results/6_bowtie2/reduntant/2_bowtie2_alignment_bam/{sample}.bam"
    conda:
        "envs/bowtie2_samtools.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["bowtie2_reduntant"]["runtime"],
        cpus_per_task = config["bowtie2_reduntant"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/6_bowtie2/reduntant/alignment_{sample}.log",
        err="log/6_bowtie2/reduntant/alignment_{sample}.err"
    shell:
        """
        mkdir -p results/6_bowtie2/reduntant/1_bam_file results/6_bowtie2/reduntant/2_bowtie2_alignment_bam
        bowtie2 --end-to-end --very-sensitive \
                -x results/6_bowtie2/reduntant/reduntant_index/ncbi_vir_reduntant_index \
                -1 {input.r1} -2 {input.r2} \
                -S {output.sam} \
                --threads {resources.cpus_per_task} > {log.out} 2> {log.err}
        
        # Samtools processing
        samtools view -@ {resources.cpus_per_task} -hbS -f 2 {output.sam} \
        | samtools sort -@ {resources.cpus_per_task} -o {output.bam} - \
            >> {log.out} 2>> {log.err}
        samtools index -@ {resources.cpus_per_task} {output.bam} \
            >> {log.out} 2>> {log.err}
        """

# Rule 10: coverm count/coverage/tpm   calculate

rule coverm_process:
    input:
        bam = "results/6_bowtie2/reduntant/2_bowtie2_alignment_bam/{sample}.bam"
    output:
        count_file = "results/6_bowtie2/reduntant/3_coverm/count/{sample}/{sample}.count.tsv",
        coverage_file = "results/6_bowtie2/reduntant/3_coverm/coverage/{sample}/{sample}.coverage.tsv",
        tpm_file = "results/6_bowtie2/reduntant/3_coverm/tpm/{sample}/{sample}.tpm.tsv"   
    conda:
        "envs/coverm.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["coverm"]["runtime"],  # minutes
        cpus_per_task = config["coverm"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/6_bowtie2/reduntant/coverm_process_{sample}.log",
        err="log/6_bowtie2/reduntant/coverm_process_{sample}.err"
    shell:
        """
        mkdir -p results/6_bowtie2/reduntant/3_coverm/count/{wildcards.sample}
        mkdir -p results/6_bowtie2/reduntant/3_coverm/coverage/{wildcards.sample}
        mkdir -p results/6_bowtie2/reduntant/3_coverm/tpm/{wildcards.sample}
        mkdir -p log/6_bowtie2/reduntant
        
        coverm contig -m count --bam-files {input.bam} \
            --min-read-percent-identity {config[coverm][min-read-percentage-identity]} \
            --min-read-aligned-percent {config[coverm][min-read-aligned-percent]} \
            --output-file {output.count_file} \
            --contig-end-exclusion {config[coverm][contig-end-exclusion]} \
            --no-zeros -t {resources.cpus_per_task} \
            >> {log.out} 2>> {log.err}
           
        coverm contig -m mean --bam-files {input.bam} \
            --min-read-percent-identity {config[coverm][min-read-percentage-identity]} \
            --min-read-aligned-percent {config[coverm][min-read-aligned-percent]} \
            --output-file {output.coverage_file} \
            --contig-end-exclusion {config[coverm][contig-end-exclusion]} \
            --no-zeros -t {resources.cpus_per_task} \
            >> {log.out} 2>> {log.err}

        coverm contig -m tpm --bam-files {input.bam} \
            --min-read-percent-identity {config[coverm][min-read-percentage-identity]} \
            --min-read-aligned-percent {config[coverm][min-read-aligned-percent]} \
            --output-file {output.tpm_file} \
            --contig-end-exclusion {config[coverm][contig-end-exclusion]} \
            --no-zeros -t {resources.cpus_per_task} \
            >> {log.out} 2>> {log.err}
            
        echo "CoverM processing completed for {wildcards.sample}."
        """
                    



