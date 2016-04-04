shell.prefix("set -euo pipefail;")
configfile: "config.yaml"


TISSUES = config["samples_pe"]
n_samples   = config["subsamplings"]["n_samples"],
sample_size = config["subsamplings"]["sample_size"],
sample_idx  = list([str(x+1) for x in range(int(n_samples[0]))])
SEED    = config["random"]["seed"] 


include: "snakefiles/folders"
include: "snakefiles/software"
include: "snakefiles/download"
include: "snakefiles/raw"
include: "snakefiles/subsampling"
include: "snakefiles/qc"
include: "snakefiles/genome"
include: "snakefiles/transcriptome"
include: "snakefiles/assembly"

# Variables



rule all:
    input:
        expand(
            gff_dna + "{tissue}_{sampling}.gff",
            tissue = config["samples_pe"],
            sampling = sample_idx
        ),
        expand(
            assembly + "{tissue}_{sampling}.fa",
            tissue = config["samples_pe"],
            sampling = sample_idx
        ),
        expand(
            blastn + "{tissue}_{sampling}.tsv",
            tissue = config["samples_pe"],
            sampling = sample_idx
        ),
        blastn + "matches.tsv",
        gff_dna + "counts.tsv"
        
            
        



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







