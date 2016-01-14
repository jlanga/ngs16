shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

# Folder variables
## QC
raw        = "data/fastq_raw"
sampled    = "data/fastq_samples"
trimmed    = "data/fastq_trimmed"
# References
references = "data/reference"
indexes    = "data/index"
## genome-based
bam_dna    = "data/bamg"
gff_dna    = "data/gffg"
gff_merged_dna = "data/gffmergedg"
## transcriptome-based
bam_rna    = "data/bamt"
gff_rna    = "data/gfft"
gff_merged_rna = "data/gffmergedt"
## de novo assembly based
assembly   = "data/assembly"
bam_ass    = "data/bama"


# Software
gzip        = config["software"]["gzip"]
trimmomatic = config["software"]["trimmomatic"]
hisat2build = config["software"]["hisat2build"]
hisat2      = config["software"]["hisat2"]
stringtie   = config["software"]["stringtie"]
bwa         = config["software"]["bwa"]
trinity     = config["software"]["trinity"]
samtools0   = config["software"]["samtools0"]

# Variables
TISSUES = config["samples_pe"]

SEED    = config["random"]["seed"] 

# Variables - Sampling
n_samples   = config["subsamplings"]["n_samples"],
sample_size = config["subsamplings"]["sample_size"],
sample_idx  = list([str(x+1) for x in range(int(n_samples[0]))])



rule all:
    input:
        expand(
            gff_merged_dna + "/{tissue}_{sampling}.gff",
            tissue = config["samples_pe"],
            sampling = sample_idx
        ),
        expand(
            gff_merged_rna + "/{tissue}_{sampling}.gff",
            tissue = config["samples_pe"],
            sampling = sample_idx
        ),
        expand(
            assembly + "/{tissue}_{sampling}.fa",
            tissue = config["samples_pe"],
            sampling = sample_idx
        )
            
        



rule clean:
    shell:
        """
        #rm -rf logs
        #rm -rf benchmarks
        #rm -rf {references} {indexes}
        rm -rf {raw} {sampled} {trimmed}
        rm -rf {bam_dna} {gff_dna} {gff_merged_dna}
        rm -rf {bam_rna} {gff_rna} {gff_merged_rna}
        """



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



rule download_gff:
    output:
        gff = references + "/annotation.gff"
    params:
        url = config["reference"]["gff"]
    log:
        "logs/download/gff.log"
    benchmark:
        "benchmarks/download/gff.json"
    shell:
        """
        wget \
            -O - \
            -o {log} \
            {params.url} |
        {gzip} -dc \
        > {output.gff}
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
        ln -s $(readlink -f {input.forward}) {output.forward} 2>  {log}
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
        {gzip} -9 > {output.reverse} ) 2>> {log}
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



rule genome_hisat2index:
    input:
        fa = references + "/dna.fa"
    output:
        files = expand(
            indexes + "/dna.{extensions}.bt2",
            extensions = "1 2 3 4 5 6 7 8".split()),
        mock  = touch(indexes + "/dna.txt")
    params:
        prefix = indexes + "/dna"
    threads:
        1
    log:
        "logs/genome/hisat2index.log"
    benchmark:
        "benchmarks/genome/hisat2index.json"
    shell:
        """
        {hisat2build} \
            {input.fa} \
            {params.prefix} \
        > {log} 2>&1 
        """



rule genome_hisat2_tissue_sampling:
    input:
        forward  = trimmed + "/{tissue}_{sampling}_1.fq.gz",
        reverse  = trimmed + "/{tissue}_{sampling}_2.fq.gz",
        unpaired = trimmed + "/{tissue}_{sampling}_u.fq.gz",
        index    = indexes + "/dna.txt"
    output:
        bam      = bam_dna + "/{tissue}_{sampling}.bam"
    params:
        index_prefix = indexes + "/dna",
        other_params = config["hisat2params"]["other"],
        sq_id        = "{tissue}",
        sq_library   = "LB:truseq_{tissue}",
        sq_platform  = "PL:Illumina",
        sq_sample    = "SM:{tissue}_{sampling}"
    threads:
        24
    log:
        "logs/genome/hisat2_{tissue}_{sampling}.log"
    benchmark:
        "benchmarks/genome/hisat2_{tissue}_{sampling}.json"
    shell:
        """
        ( {hisat2} \
            --rg-id {params.sq_id} \
            --rg {params.sq_library} \
            --rg {params.sq_platform} \
            --rg {params.sq_sample} \
            --threads {threads} \
            {params.other_params} \
            -x {params.index_prefix} \
            -1 {input.forward} \
            -2 {input.reverse} \
            -U {input.unpaired} |
        {samtools0} view \
            -S \
            -h \
            -u \
            - |
        {samtools0} rmdup \
            -s - - |
        {samtools0} sort \
            -@ {threads} \
            -l 9 \
            -o \
            - \
            - \
        > {output.bam} ) \
        2>  {log}
        """



rule genome_stringtie_tissue_sampling:
    input:
        bam = bam_dna + "/{tissue}_{sampling}.bam",
        gff = references + "/annotation.gff"
    output:
        gff = gff_dna + "/{tissue}_{sampling}.gff"
    threads:
        1
    params:
        label = "{tissue}_p{sampling}"
    log:
        "logs/genome/stringtie_{tissue}_{sampling}.log"
    benchmark:
        "benchmarks/genome/stringtie_{tissue}_{sampling}.json"
    shell:
        """
        {stringtie} \
            {input.bam} \
            -G {input.gff} \
            -l {params.label} \
            -o {output.gff} \
            -p {threads} \
        2> {log}
        """



def genome_merge_gff_files(wildcards):
    partitions = list([str(x + 1) for x in range(int(wildcards.sampling))])
    gffs = [gff_dna + "/" + wildcards.tissue + "_" + partition +  ".gff" for partition in partitions]
    return gffs 


rule genome_stringtie_merge_tissue_sampling:
    input:
        samples_gff = genome_merge_gff_files,
        reference_gff = references + "/annotation.gff"
    output:
        gff = gff_merged_dna + "/{tissue}_{sampling}.gff"
    threads:
        1
    params:
        label = "{tissue}_{sampling}"
    log:
        "logs/genome/stringtie_merge_{tissue}_{sampling}.log"
    benchmark:
        "benchmarks/genome/stringtie_merge_{tissue}_{sampling}.json"
    shell:
        """
        {stringtie} --merge \
            -G {input.reference_gff} \
            -o {output.gff} \
            -l {params.label} \
            {input.samples_gff} \
        2> {log}
        """
    
    

rule transcriptome_hisat2index:
    input:
        fa = references + "/cdna.fa"
    output:
        files = expand(
            indexes + "/cdna.{extensions}.bt2",
            extensions = "1 2 3 4 5 6 7 8".split()),
        mock  = touch(indexes + "/cdna.txt")
    params:
        prefix = indexes + "/cdna"
    threads:
        1
    log:
        "logs/transcriptome/hisat2index.log"
    benchmark:
        "benchmarks/transcriptome/hisat2index.json"
    shell:
        """
        {hisat2build} \
            {input.fa} \
            {params.prefix} \
        > {log} 2>&1 
        """



rule transcriptome_hisat2t_tissue_sampling:
    input:
        forward  = trimmed + "/{tissue}_{sampling}_1.fq.gz",
        reverse  = trimmed + "/{tissue}_{sampling}_2.fq.gz",
        unpaired = trimmed + "/{tissue}_{sampling}_u.fq.gz",
        index    = indexes + "/cdna.txt"
    output:
        bam      = bam_rna + "/{tissue}_{sampling}.bam"
    params:
        index_prefix = indexes + "/cdna",
        other_params = config["hisat2params"]["other"],
        sq_id        = "{tissue}",
        sq_library   = "LB:truseq_{tissue}",
        sq_platform  = "PL:Illumina",
        sq_sample    = "SM:{tissue}_{sampling}"
    threads:
        24
    log:
        "logs/transcriptome/hisat2_{tissue}_{sampling}.log"
    benchmark:
        "benchmarks/transcriptome/hisat2_{tissue}_{sampling}.json"
    shell:
        """
        ( {hisat2} \
            --rg-id {params.sq_id} \
            --rg {params.sq_library} \
            --rg {params.sq_platform} \
            --rg {params.sq_sample} \
            --threads {threads} \
            {params.other_params} \
            -x {params.index_prefix} \
            -1 {input.forward} \
            -2 {input.reverse} \
            -U {input.unpaired} |
        {samtools0} view \
            -S \
            -h \
            -u \
            - |
        {samtools0} rmdup \
            -s - - |
        {samtools0} sort \
            -@ {threads} \
            -l 9 \
            -o \
            - \
            - \
        > {output.bam} ) \
        2>  {log}
        """



rule transcriptome_stringtie_tissue_sampling:
    input:
        bam = bam_rna + "/{tissue}_{sampling}.bam",
    output:
        gff = gff_rna + "/{tissue}_{sampling}.gff"
    threads:
        1
    params:
        label = "{tissue}_p{sampling}"
    log:
        "logs/transcriptome/stringtie_{tissue}_{sampling}.log"
    benchmark:
        "benchmarks/transcriptome/stringtie_{tissue}_{sampling}.json"
    shell:
        """
        {stringtie} \
            {input.bam} \
            -l {params.label} \
            -o {output.gff} \
            -p {threads} \
        2> {log}
        """



def transcriptome_merge_gff_files(wildcards):
    partitions = list([str(x + 1) for x in range(int(wildcards.sampling))])
    gffs = [gff_rna + "/" + wildcards.tissue + "_" + partition +  ".gff" for partition in partitions]
    return gffs 



rule transcriptome_stringtie_merge_tissue_sampling:
    input:
        samples_gff = transcriptome_merge_gff_files
    output:
        gff = gff_merged_rna + "/{tissue}_{sampling}.gff"
    threads:
        1
    params:
        label = "{tissue}_{sampling}"
    log:
        "logs/transcriptome/stringtie_merge_{tissue}_{sampling}.log"
    benchmark:
        "benchmarks/transcriptome/stringtie_merge_{tissue}_{sampling}.json"
    shell:
        """
        {stringtie} --merge \
            -o {output.gff} \
            -l {params.label} \
            {input.samples_gff} \
        2> {log}
        """



def list_fastqs_to_assemble_forward(wildcards):
    partitions = list([str(x + 1) for x in range(int(wildcards.sampling))])
    forwards = [
        trimmed + "/" + wildcards.tissue + "_" + partition + "_1.fq.gz"
        for partition in partitions
    ]
    return forwards

def string_fastqs_to_assemble_forward(wildcards):
    return ",".join(list(list_fastqs_to_assemble_forward(wildcards)))

def list_fastqs_to_assemble_reverse(wildcards):
    partitions = list([str(x + 1) for x in range(int(wildcards.sampling))])
    reverses = [
        trimmed + "/" + wildcards.tissue + "_" + partition + "_2.fq.gz"
        for partition in partitions
    ]
    return reverses

def string_fastqs_to_assemble_reverse(wildcards):
    return ",".join(list(list_fastqs_to_assemble_reverse(wildcards)))


rule assembly_run_trinity_tissue_sampling:
    input:
        left  = list_fastqs_to_assemble_forward,
        right = list_fastqs_to_assemble_reverse
    output:
        fasta = assembly + "/{tissue}_{sampling}.fa"
    threads:
        24
    params:
        left  = string_fastqs_to_assemble_forward,
        right = string_fastqs_to_assemble_reverse, 
        memory= config["trinity_params"]["memory"],
        outdir= assembly + "/trinity_{tissue}_{sampling}"
    log:
        "logs/assembly/{tissue}_{sampling}.log"
    benchmark:
        "benchmarks/assembly/{tissue}_{sampling}.log"
    shell:
        """
        ./bin/Trinity \
            --seqType fq \
            --max_memory {params.memory} \
            --left {params.left} \
            --right {params.right} \
            --CPU {threads} \
            --full_cleanup \
            --output {params.outdir} \
        > {log}
        
        mv {params.outdir}.Trinity.fasta {output.fasta}
        
        #rm {input.left}.readcount {input.right}.readcount
        """



