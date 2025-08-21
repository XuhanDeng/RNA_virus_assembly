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
        
        # Step 6: Custom RDRP database and analysis
        "results/5_blastx/database/diamond/custom_rdrp.dmnd",
        "results/7_hmmer/database/custom_modules.hmm",
        
        # Step 7: BLASTX analysis  
        expand("results/5_blastx/blastx/{sample}/{sample}_blastx.tbl", sample=config["samples"]),
        
        # Step 8: HMM search results
        "results/7_hmmer/hmmer_hits.tbl"
        





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
        fastp --thread {config[fastp][threads]} \
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
        ribodetector_cpu -t {config[ribodetector][threads]} \
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
                  -t {config[spades][threads]} -m {resources.mem_mb_per_cpu} \
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

# rule 5 established viral RDRP protein database
# NCBI ref viral database downlaed on 2025-08-15(3,105 sequences)
# well-currated rdrp database from "Using artificial intelligence to document the hidden RNA virosphere" 
# RDRP database from "Using artificial intelligence to document the hidden RNA virosphere"  identification results (513,134 sequences)
#cd-hit was used for a threshold with 90% identity 

rule create_rdrp_protein_database:
    input:
        cell_well_currated = "database/rdrp_all/Cell_well/5979_known_viral_RdRPs.fasta",
        cell_identified = "database/rdrp_all/Cell_identified/cell_identified.faa",
        NCBI_ref = "database/rdrp_all/ncbi_0815/NCBI_ref_rdrp_0815.faa"
    output:
        db = "results/5_blastx/database/diamond/custom_rdrp.dmnd",
        clustered_fasta = "results/5_blastx/database/90/custom_rdrp.90.fasta"
    conda:
        "envs/diamond-seqkit-cdhit.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],
        runtime = 6000,  # minutes (12 hours)
        cpus_per_task = 32,
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/create_rdrp_protein_db/makedb.log",
        err="log/create_rdrp_protein_db/makedb.err"
    shell:
        """
        mkdir -p log/create_rdrp_protein_db results/5_blastx/database/90 results/5_blastx/database/diamond
        cat {input.cell_well_currated} {input.cell_identified} {input.NCBI_ref} > results/5_blastx/database/merge_custom_rdrp.fasta
        
        cd-hit -i results/5_blastx/database/merge_custom_rdrp.fasta -o {output.clustered_fasta} -c 0.9 -n 5 -g 1 -G 0 -aS 0.8 -d 0 -p 1 -T {resources.cpus_per_task} -M 0 -B 1
        
        diamond makedb --in {output.clustered_fasta} --db {output.db} > {log.out} 2> {log.err}
        """




# rule6 establish HMM database for RDRP
# HMM database is built from the clustered RDRP protein sequences with 40% identity
rule hmmdatabase_establishment:
    input:
        clustered_fasta = "results/5_blastx/database/90/custom_rdrp.90.fasta"
    output:
        clustered_fasta_40 = "results/7_hmmer/database/cdhit/40/cdhit40.fasta",
        hmm_profile = "results/7_hmmer/database/custom_modules.hmm"
    conda: "envs/cdhit-biopython-mafft-hmmer.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],
        runtime = 6000,  # minutes
        cpus_per_task = 16,
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/7_hmmer/hmmdatabase_establishment.log",
        err="log/7_hmmer/hmmdatabase_establishment.err"
    shell:
        """
        mkdir -p results/7_hmmer/database/cdhit/40 log/7_hmmer
        mkdir -p results/7_hmmer/database/{{fasta,mafft,hmm}}
        
        # First clustering at 40% identity using PSI-CD-HIT
        perl {config[software][psi-cd-hit]} -i {input.clustered_fasta} -o {output.clustered_fasta_40} -c 0.4 -ce 1e-3 -aS 0.5 -G 1 -g 1 -exec local -para 16 -blp 16
        
        # Parse clusters and extract sequences for each cluster
        python {config[software][ParseCDHIT]} results/7_hmmer/database/cdhit/40/cdhit40.fasta.clstr {input.clustered_fasta}
        mv *.fasta results/7_hmmer/database/fasta/
        
        # Align each cluster and build HMM profiles
        for file in results/7_hmmer/database/fasta/*.fasta; do
            if [ -f "$file" ]; then
                basename=$(basename "$file" .fasta)
                # Align sequences in cluster
                mafft --thread {resources.cpus_per_task} --localpair --maxiterate 1000 "$file" > results/7_hmmer/database/mafft/"$basename".aln.faa
                # Build HMM profile for cluster
                hmmbuild results/7_hmmer/database/hmm/"$basename".hmm results/7_hmmer/database/mafft/"$basename".aln.faa
            fi
        done
        
        # Combine all HMM profiles
        cat results/7_hmmer/database/hmm/*.hmm > {output.hmm_profile}
        hmmpress {output.hmm_profile}
        """



# Rule 7: Run BLASTX to find RDRP sequences in assemblies
rule blastx_rdrp:
    input:
        fasta = "results/4_reformat/{sample}_reformat.fasta",
        db = "results/5_blastx/database/diamond/custom_rdrp.dmnd"
    output:
        tbl = "results/5_blastx/blastx/{sample}/{sample}_blastx.tbl",
        unique_ids = "results/5_blastx/blastx/{sample}/{sample}_blastx.unique_ids.txt"
    conda:
        "envs/diamond-seqkit-cdhit.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],
        runtime = 6000,  # minutes
        cpus_per_task = config["blastx"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/5_blastx/blastx/{sample}.log",
        err="log/5_blastx/blastx/{sample}.err"
    shell:
        """
        mkdir -p results/5_blastx/blastx/{wildcards.sample} log/5_blastx/blastx
        
        # Filter sequences with minimum length of 600 bp
        seqkit seq -m 600 {input.fasta} > results/5_blastx/blastx/{wildcards.sample}/{wildcards.sample}_filtered.fasta
        
        # Run DIAMOND BLASTX
        diamond blastx --query results/5_blastx/blastx/{wildcards.sample}/{wildcards.sample}_filtered.fasta \
            --db {input.db} \
            --out {output.tbl} \
            --outfmt 6 qseqid sseqid qlen length pident evalue \
            --evalue 1E-5 \
            --threads {config[blastx][threads]} \
            --max-target-seqs 1 \
            --very-sensitive \
            --quiet > {log.out} 2> {log.err}
        
        # Extract unique sequence IDs from BLASTX hits
        cat {output.tbl} | cut -f1 | sort -u > {output.unique_ids}
        """


#Rule 8 : orfipy to get ORFs for each sample
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
        orfipy {input.fasta} --dna {output.orfs} --min 10 --max 10000 --procs {config[orfipy][threads]} --outdir results/5_orfinder/{wildcards.sample}/ > {log.out} 2> {log.err}
        """




rule hmmsearch_orfs:
    input:
        hmmdb = "results/7_hmmer/database/custom_modules.hmm",
        orfs = "results/5_orfinder/merged_orfs.fasta"
    output:
        tbl = "results/7_hmmer/hmmer_hits.tbl",
        domtbl = "results/7_hmmer/hmmer_dom.tbl"
    conda: "envs/cdhit-biopython-mafft-hmmer.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],
        runtime = 180,  # minutes
        cpus_per_task = 16,
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/7_hmmer/hmmsearch.log",
        err="log/7_hmmer/hmmsearch.err"
    shell:
        """
        mkdir -p results/7_hmmer log/7_hmmer
        hmmscan --cpu {resources.cpus_per_task} --tblout {output.tbl} --domtblout {output.domtbl} {input.hmmdb} {input.orfs} > {log.out} 2> {log.err}
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