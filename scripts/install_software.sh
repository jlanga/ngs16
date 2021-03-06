#!/usr/bin/env bash

mkdir -p src
pushd src

# hisat2-2.0.1
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.1-beta-Linux_x86_64.zip
unzip hisat2-2.0.1-beta-Linux_x86_64.zip
pushd hisat2-2.0.1-beta
cp extract_* hisat* simulate_reads.py ../../bin/
popd

# samtools v0.1.20
wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar -xvf samtools-0.1.19.tar.bz2
pushd samtools-0.1.19
make -j 8
cp samtools ../../bin/
popd

# bwa-0.7.12
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2
tar -xvf bwa-0.7.12.tar.bz2
pushd bwa-0.7.12
make -j 8
cp bwa ../../bin/

# trimmomatic-0.35
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip
unzip Trimmomatic-0.35.zip
cp Trimmomatic-0.25/trimmomatic-0.35.jar ../bin/

# stringtie-1.2.0
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.2.0.Linux_x86_64.tar.gz
tar xvf stringtie-1.2.0.Linux_x86_64.tar.gz
cp stringtie-1.2.0.Linux_x86_64/stringtie ../bin/stringtie

# trinity-2.1.1
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz
tar xvf v2.1.1.tar.gz
pushd trinityrnaseq-2.1.1
make -j 8 
make -j 8 plugins
make test
popd
ln -s ../src/trinityrnaseq-2.1.1/Trinity ../bin/Trinity

popd
