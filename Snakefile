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

# Variables
TISSUES = config["samples_pe"]

SEED = config["random"]["seed"] 

# Variables - Sampling
n_samples   = config["subsamplings"]["n_samples"],
sample_size = config["subsamplings"]["sample_size"],
sample_idx  = list([str(x+1) for x in range(int(n_samples[0]))])


rule all:
    input:
        references + "/cdna.fa",
        references + "/dna.fa",
        expand(
            trimmed + "/{tissue}_{sampling}_{pair}.fq.gz",
            tissue = config["samples_pe"],
            sampling = sample_idx,
            pair = "1 2 u".split()
        )


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



rule raw_make_links_pe:
    input:
        forward= lambda wildcards: config["samples_pe"][wildcards.tissue]["forward"],
        reverse= lambda wildcards: config["samples_pe"][wildcards.tissue]["reverse"]
    output:
        forward= raw + "/{tissue}_1.fq.gz",
        reverse= raw + "/{tissue}_2.fq.gz"
    threads:
        1
    log:
        "logs/raw/make_links_pe_{tissue}.log"
    benchmark:
        "benchmarks/raw/make_links_pe_{tissue}.json"
    shell:
        """
        ln -s $(readlink -f {input.forward}) {output.forward} 2> {log}
        ln -s $(readlink -f {input.reverse}) {output.reverse} 2>> {log}
        """


rule subsampling_subsample_integers_sample:
    output:
        expand(
            sampled + "/{{tissue}}_{sampling}.idx",
            sampling = sample_idx
        )
    threads:
        1
    params:
        prefix = sampled + "/{tissue}",
        set_size = lambda wildcards: config["samples_pe"][wildcards.tissue]["reads"],
        n_samples   = n_samples,
        sample_size = sample_size,
        seed = SEED,
    log:
        "logs/subsampling/subsample_integers_{tissue}.log"
    benchmark:
        "benchmarks/subsampling/subsample_integers_{tissue}.json"
    shell:
        """
        python3 scripts/subsample_integers.py \
            {params.prefix} \
            {params.set_size} \
            {params.n_samples} \
            {params.sample_size} \
            {params.seed} \
        2> {log}
        """



rule subsampling_filter_fastq_by_record_number_sample:
    input:
        forward = raw + "/{tissue}_1.fq.gz",
        reverse = raw + "/{tissue}_2.fq.gz",
        index   = sampled + "/{tissue}_{sampling}.idx"
    output:
        forward = sampled + "/{tissue}_{sampling}_1.fq.gz",
        reverse = sampled + "/{tissue}_{sampling}_2.fq.gz"
    threads:
        2
    params:
    
    log:
        "logs/subsampling/filter_fastq_by_record_number_{tissue}.log"
    benchmark:
        "benchmarks/subsampling/filter_fastq_by_record_number_{tissue}.json"
    shell:
        """
        ( {gzip} -dc {input.forward} |
        python3 scripts/filter_fastq_by_record_number.py \
            {input.index} |
        {gzip} -9 > {output.forward} ) 2> {log}
        
        ( {gzip} -dc {input.reverse} |
        python3 scripts/filter_fastq_by_record_number.py \
            {input.index} |
        {gzip} -9 > {output.reverse} ) 2> {log}
        """



rule qc_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and 
    remove low quality regions and reads.
    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
    Outputs _3 and _4 are compressed with the builtin compressor from 
    Trimmomatic. Further on they are catted and compressed with gzip/pigz 
    (level 9).
    Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ
    header. It is done posterior to the trimming since the output comes 
    slower than the input is read.
    Number of threads used:
        4 for trimmomatic
        2 for gzip inputs
        2 for gzip outputs
        Total: 8
    """
    input:
        forward    = sampled + "/{tissue}_{sampling}_1.fq.gz",
        reverse    = sampled + "/{tissue}_{sampling}_2.fq.gz"
    output:
        forward    = trimmed + "/{tissue}_{sampling}_1.fq.gz",
        reverse    = trimmed + "/{tissue}_{sampling}_2.fq.gz",
        unpaired   = trimmed + "/{tissue}_{sampling}_u.fq.gz"
    params:
        unpaired_1 = trimmed + "/{tissue}_{sampling}_3.fq.gz",
        unpaired_2 = trimmed + "/{tissue}_{sampling}_4.fq.gz",
        adaptor    = lambda wildcards: config["samples_pe"][wildcards.tissue]["adaptor"],
        phred      = lambda wildcards: config["samples_pe"][wildcards.tissue]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    benchmark:
        "benchmarks/qc/trimmomatic_pe_{tissue}.json"
    log:
        "logs/qc/trimmomatic_pe_{tissue}.log" 
    threads:
        8
    shell:
        """
        {trimmomatic} PE \
            -threads {threads} \
            -{params.phred} \
            <({gzip} -dc {input.forward} ) \
            <({gzip} -dc {input.reverse} ) \
            >(cut -f 1 -d " " | {gzip} -9 > {output.forward} ) \
            {params.unpaired_1} \
            >(cut -f 1 -d " " | {gzip} -9 > {output.reverse} ) \
            {params.unpaired_2} \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log}
            
        zcat {params.unpaired_1} {params.unpaired_2} |
        cut -f 1 -d " " |
        {gzip} -9 > {output.unpaired}
        
        rm {params.unpaired_1} {params.unpaired_2}
        """
