shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

# Folder variables
raw        = "data/fastq_raw"
sampled    = "data/fastq_samples"
trimmed    = "data/fastq_trimmed"
bam_dna    = "data/bamg"
bam_rna    = "data/bamt"
assembly   = "data/assembly"
bam_ass    = "data/bama"
references = "data/reference"
indexes    = "data/index"

# Software
gzip        = config["software"]["gzip"]
trimmomatic = config["software"]["trimmomatic"]
hisat2      = config["software"]["hisat2"]
bwa         = config["software"]["bwa"]
trinity     = config["software"]["trinity"]
samtools0   = config["software"]["samtools0"]




rule all:
    


rule download_dna:
    output:
        fa = references + "/dna.fa"
    threads:
        1
    params:
        url = config["reference"]["dna"]
    log:
        "logs/download/dna.log"
    benchmark:
        "benchmarks/download/dna.json"
    shell:
        """
        wget                \
            -O  -           \
            -o  {log}       \
            {params.url}    |
        {gzip} -dc          \
        > {output.fa}
        """  



rule download_cdna:
    output:
        fa = references + "/cdna.fa"
    threads:
        1
    params:
        url = config["reference"]["cdna"]
    log:
        "logs/download/cdna.log"
    benchmark:
        "benchmarks/download/cdna.json"
    shell:
        """
        wget                \
            -O  -           \
            -o  {log}       \
            {params.url}    |
        {gzip} -dc          \
        > {output.fa}
        """
