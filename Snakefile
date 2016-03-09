shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

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
TISSUES = config["samples_pe"]







rule all:
    input:
        #expand(
        #    gff_dna + "/{tissue}_{sampling}.gff",
        #    tissue = config["samples_pe"],
        #    sampling = sample_idx
        #),
        #expand(
        #    gff_merged_rna + "/{tissue}_{sampling}.gff",
        #    tissue = config["samples_pe"],
        #    sampling = sample_idx
        #),
        expand(
            assembly + "{tissue}_{sampling}.fa",
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







