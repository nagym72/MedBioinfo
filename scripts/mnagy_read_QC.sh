#!/bin/bash

echo "script start: download and initial sequencing read quality control"
date

#fetch out corresponding patients that belong to mnagy
#then direct output to .txt file

sqlite3 --noheader --csv /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db \
"select run_accession from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='mnagy';" > mnagy_run_accessions.txt

#load fastq module
module load fastq

# for each run pipe into fastq-dump with xargs -I to take multiple args
# taken from the .txt file. These contain the Identifiers required
# to download corresponding fastq files.
# finally gzip all outputfiles with --gzip and specify outdir with 
# -- outdir
# no multithreading here.
# --split-3 in order to generate fastq files for :

#forward, reverse and in case of no direction undirected outfiles.
#load sra-tools

module load sra-tools

#make dir
mkdir ../data/sra_fastq

cat mnagy_run_accessions.txt | srun --cpus-per-task=1 \
--time=00:30:00 xargs \
fastq-dump -I --gzip --outdir ../data/sra_fastq/ \
--disable-multithreading --split-3

#load seqkit module

#seqkit load module
module load seqkit

#seqkit to get stats, this numbers coincides 
# with zcat <file> | grep "^@" | wc -l

#check how many reads present.
srun --cpus-per-task=1 --time=00:05:00 seqkit stats \
../data/sra_fastq/*.fastq.gz


# check for duplicates in forward runs.
# seqkit rmdup -s to check through seqs
# -D to direct seq ids to a specified outfile.

cat mnagy_run_accessions.txt | xargs -I {} zcat ../data/sra_fastq/{}_1.fastq.gz |\
srun --cpus-per-task=1 --time=00:05:00 seqkit rmdup -s -D \
dupe_seq_ids_forward_runs.txt 

#check for duplicates in reverse runs.
cat mnagy_run_accessions.txt | xargs -I {} zcat ../data/sra_fastq/{}_2.fastq.gz | \
srun --cpus-per-task=1 --time=00:05:00 seqkit rmdup -s -D \
dupe_seq_ids_reverse_runs.txt


#search for adapters. full were not found.
#so lets search for shorter sequences.

#again for all forward and reverse runs.


# -i ignore case 
# -d pattern contains degenerate motif
# -p searchpattern AGATCGGAAGAGC

cat mnagy_run_accessions.txt | xargs -I {} zcat ../data/sra_fastq/{}_1.fastq.gz | \
srun --cpus-per-task=1 --time=00:05:00 seqkit locate -i -d -p \
AGATCGGAAGAGC -o adapter_search.txt

                                                    
#fastqc quality run.

mkdir ../analyses/fastqc

srun --cpus-per-task=4 --time=00:15:00 xargs -I{} -a mnagy_run_accessions.txt fastqc
 --o ../analyses/fastqc --noextract --threads 4 ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz


#trimmomatic additional work.

#check if adapters were trimmed.
#expected result:already trimmed.

srun --cpus-per-task=4 --time=00:10:00 trimmomatic PE \
../data/sra_fastq/ERR6913137_1.fastq.gz \
../data/sra_fastq/ERR6913137_2.fastq.gz \
../data/sra_fastq/trimmomatic_test/ERR6913137_1_paired_trimmed.fastq.gz \
../data/sra_fastq/trimmomatic_test/ERR6913137_1_unpaired_trimmed.fastq.gz \
../data/sra_fastq/trimmomatic_test/ERR6913137_2_paired_trimmed.fastq.gz \
../data/sra_fastq/trimmomatic_test/ERR6913137_2_unpaired_trimmed.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > trimmomatic_ERR6913137_merged.txt

#outout is 69M 71M for paired end trimmed. This is exactly the same for the inputs and statistic told us
#that we only trimmed <1% as expected because it was already trimmed.

#load flash2
module load flash2


#run flash2 on all runs.

srun --cpus-per-task=2 --time=00:30:00 xargs -a mnagy_run_accessions.txt \
-n 1 -I{} flash2 --threads=2 -z \
--output-prefix={}.flash --output-directory=../data/merged_pairs/ \
../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 2>&1 \
| tee -a mnagy_flash2.log


#check for contamination of PhiX phage.
#after setting up efetch through wget.
#install efetch
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

#store fasta file
efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna

#check if the file is correct.
#head -n 5 ../data/reference_seqs/PhiX_NC_001422.fna

#load bowtie2
module load bowtie2

#make bowtie2 index db dir
mkdir ../data/bowtie2_DBs

#make index db
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_DBs/PhiX_bowtie2_DB

#to see what files are created.
#ls ../data/bowtie2_DBs/

#new dir for analyses.
mkdir ../analyses/bowtie

#check for phage contamination.
#bowtie2 -x specify DB index files (base name is sufficient to fetch
#all index files .1 .2 etc
# -S for SAM outformat.
#threads matching cpus-per-task
# no-unal : suppress unaligned sam outfile records.
# 2>&1 collect stdout + stderr

srun --cpus-per-task=8 --time=00:10:00 bowtie2 -x ../data/bowtie2_DBs/PhiX_bowtie2_DB \
-U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz \
-S ../analyses/bowtie/mnagy_botwie_merged2PhiX.sam --threads=8 \
--no-unal 2>&1 | tee analyses/bowtie/mnagy_botwie_merged2PhiX.log


#now with SarsCov2

efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_045512.fna

#build index DB for SC2

srun bowtie2-build -f ../data/reference_seqs/SC2_045512.fna \
..data/bowtie2_DBs/SC2_bowtie2_DB

#lets run bowtie2

srun --cpus-per-task=8 --time=00:03:00 bowtie2 -x ../data/bowtie2_DBs/ \
-U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S \
../analyses/mnagy_merged2SC2.sam --threads=8 \
--no-unal 2>&1 | tee ../analyses/bowtie/mnagy_botwie_merged2SC2.log

#load multiqc module
module load multiqc

#combine muiltiple reports from fastqc/flash2/bowtie
srun multiqc --force --title "michael nagy sample sub-set" \
--comment "Summary of Assignment 10 for medbioinfo" \
../data/merged_pairs/ ../analyses/fastqc/ ./mnagy_flash2.log ../analyses/bowtie/


#on my local machine to transfer multiqc outfile summary html file.

#this will transfer html files to the current dir I am on my machine ( same as previous scp destination)
#scp mnagy@core.cluster.france-bioinformatique.fr:/shared/home/mnagy/medbioinfo_folder/michael/MedBioinfo/scripts/*.html \
#./ 

#ls on my local:
#ERR6913137_1_fastqc.html  ERR6913233_2_fastqc.html
#ERR6913137_2_fastqc.html  ERR6913255_1_fastqc.html
#ERR6913159_1_fastqc.html  ERR6913255_2_fastqc.html
#ERR6913159_2_fastqc.html  ERR6913262_1_fastqc.html
#ERR6913166_1_fastqc.html  ERR6913262_2_fastqc.html
#ERR6913166_2_fastqc.html  ERR6913269_1_fastqc.html
#ERR6913173_1_fastqc.html  ERR6913269_2_fastqc.html
#ERR6913173_2_fastqc.html  michael-nagy-sample-sub-set_multiqc_report.html
#ERR6913233_1_fastqc.html

# UPON INSPECTION OF THE HTML FILE:

# IT STATES THERE IS NO PHIX contamination (0.0%) which is good.
# IT ALSO STATES THERE IS NO COVSARS (0.0%) even though it states 
# SE MAPPED UNIQUELY: 419 (so not 0.0% but it seems to round).
# A BIT MISSLEADING.

# FASTQC OUTPUT:
# ERR6913137_1 + 2 had 7% duplicated read counts.
# ERR6913233_1 + 2 had 58.4% duplicated read counts.
# ERR6913262_1 + 2 had 50% duplicated read counts.
# #$$6913269_1 + 2 had 41.9% duplicated read counts.

#VERY USEFUL TOOL and the html file is pretty nice to inspect.

# I made the mistake to only inspect 1 patient (both runs merged) and concluded
# that none of my samples contained duplicates (which is wrong.)
# Here I see that this 1 patient sample is full of duplicates.
# Also the others have duplicates (albeit <10%)

#BASEQUALITY

# shows what I have already seen at seqkit stats that low quality bases were trimmed.

#PER SEQUENCE GC CONTENT

#seemed to have failed in 10 cases (?why?)
#Sequence duplication level confirmed what I have seen above regarding duplicates.
 
#NO adapter found (so they were trimmed), this was seen before also nicely.

#also nice from them to include a 6min youtube video to look https://www.youtube.com/watch?v=qPbIlO_KWN0

