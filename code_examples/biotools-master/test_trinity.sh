#!/bin/sh

# Test Trinity and RSEM

mkdir -p ~/src
cd ~/src
git clone 'https://github.com/bli25wisc/RSEM.git'
cd RSEM
make
make ebseq
export PATH=$PATH:~/src/RSEM

cd ..
git clone 'https://github.com/trinityrnaseq/trinityrnaseq.git'
cd trinityrnaseq/
make clean
make

cd sample_data/test_Trinity_Assembly
./runMe.sh 
