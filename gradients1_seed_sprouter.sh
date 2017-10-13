#!/bin/bash

# gradients1_seed_sprouter.sh

# usage: seed_sprouter.sh myfavoritegene

# Run this script in a fresh directory with a single seed sequence for myfavoritegene.
# Script will look for myfavoritegene.fasta and run analysis from there
# The header for myfavoritegene.fasta should be 'friendly' with a taxid included. Example:
# >Fragilariopsis_cylindrus_CCMP1102_gbOEU23745_tax635003

# USAGE: seed_sprouder.sh [gene_id] [NCORES] [BITSCORE_CUTOFF] [USEARCH_THRESH] [PENDANT_LENGTH_CUTOFF]
# EXAMPLE: seed_sprouder.sh ftn_v2 16 45 95 0.7
# these two are manual inputs:
GENE=$1
NCORES=$2 # e.g., 16 for gross and 8 for match
BITSCORE_CUTOFF=$3 # original default value was 35
USEARCH_THRESH=$4 # original default was 97
PENDANT_LENGTH_CUTOFF=$5 # original default was 0.7

BLOOM_DIR="/mnt/nfs/ryan/Gradients1/runs/$GENE"
JGI_PATH="/mnt/nfs/ryan/JGI/JGI_microeuks.aa.fasta"
NCBI_PATH="/mnt/nfs/ryan/NCBI/extra_opisthokonts/marine_opisthokonts.fasta"
MARINEREF_PATH="/mnt/nfs/ryan/MarineRef/MarineRefII_ALL_SEQS.fa"
MARINEREF_INFO="/mnt/nfs/ryan/MarineRef/MarineRefII_seqIDinfo.csv"
USEARCH_PATH="/mnt/nfs/home/rgrous83/bin/usearch"
PPLACER_PATH="/mnt/nfs/home/rgrous83/bin/pplacer"
SUBJECT_DIR="/mnt/nfs/ryan/Gradients1/mORFeus_v2"
SUBJECT_LIST="/mnt/nfs/ryan/Gradients1/mORFeus_v2/gradients1_morfeus_handles.txt"
TAX_DB="/mnt/nfs/home/rgrous83/NCBI/taxonomy.db"

cd $BLOOM_DIR
mkdir jackhmmer; cd jackhmmer

# jackhmmer to recruit homologs from MarRef
# this should be fixed at <10 cores so we're leaving it hard-coded at 8.
jackhmmer -T $BITSCORE_CUTOFF --incT $BITSCORE_CUTOFF --cpu 8 -o $GENE.jackhmmer.log -A $GENE.jackhmmer.sto --tblout $GENE.jackhmmer.tab ../$GENE.fasta $MARINEREF_PATH

# switch to fasta and back to change filename
seqmagick convert $GENE.jackhmmer.sto $GENE.jackhmmer.fasta

# remove the seed sequence from the result aln:
head -1 ../$GENE.fasta | sed 's/^>//g' | awk -F" " {'print $1'} > remove_me.txt
seqmagick mogrify --exclude-from-file remove_me.txt $GENE.jackhmmer.fasta
rm remove_me.txt

# remove duplicate matches
/mnt/nfs/ryan/scripts/dedupe_fasta_by_defline.py $GENE.jackhmmer.fasta
# rename the results
python /mnt/nfs/home/rgrous83/scripts/Rename_MarineRefII_Seqs_fixed.py -I $MARINEREF_INFO -S $GENE.jackhmmer.dd.fasta

# build the hmm profile before we lose the aligment in the ucluster step:
hmmbuild $GENE.hmm $GENE.jackhmmer.fasta

# recruit sequences from diel1
mkdir ../hmmer_env; cp $GENE.hmm ../hmmer_env; cd ../hmmer_env

function hmmer_time {
SUBJECT_FASTA="$1".6tr.orfs40.fasta.gz
hmmsearch -T $BITSCORE_CUTOFF --incT $BITSCORE_CUTOFF --cpu $NCORES --tblout $GENE.$1.hmm_out.tab -A $GENE.$1.query.sto $GENE.hmm $SUBJECT_DIR/$SUBJECT_FASTA
}

for subject in $(cat $SUBJECT_LIST); do
hmmer_time $subject
done

# now we'll go and build our tree + refpkg
cd $BLOOM_DIR
mkdir hmmer_ref; cp jackhmmer/named_$GENE.jackhmmer.dd.fasta hmmer_ref; cd hmmer_ref

# use the MarRef profile to recruit from NCBI and JGI to build the tree:
hmmsearch -T $BITSCORE_CUTOFF --incT $BITSCORE_CUTOFF --cpu $NCORES --tblout $GENE.JGI.hmm_out.tab -A "$GENE".JGI.query.sto ../jackhmmer/$GENE.hmm $JGI_PATH
hmmsearch -T $BITSCORE_CUTOFF --incT $BITSCORE_CUTOFF --cpu $NCORES --tblout $GENE.NCBI.hmm_out.tab -A "$GENE".NCBI.query.sto ../jackhmmer/$GENE.hmm $NCBI_PATH

# convert the other recruits to FASTA for naming:
for source in JGI NCBI; do
	seqmagick convert $GENE.$source.query.sto $GENE.$source.query.fasta
	$USEARCH_PATH/usearch -sortbylength $GENE.$source.query.fasta -fastaout $GENE.$source.query.sorted.fasta
done
$USEARCH_PATH/usearch -sortbylength named_$GENE.jackhmmer.dd.fasta -fastaout $GENE.jackhmmer.sorted.fasta

fnames2.py -tmc $GENE.JGI.query.sorted.fasta
/mnt/nfs/ryan/scripts/dedupe_fasta_by_defline.py $GENE.JGI.query.sorted.fn.fasta
fnames2.py -tmc $GENE.NCBI.query.fasta
/mnt/nfs/ryan/scripts/dedupe_fasta_by_defline.py $GENE.NCBI.query.fn.fasta

# note that usearch drops the description so we lose tax_ids in the NCBI samples; so we're not using the sorted ones.

# Add up all the reference sequences with JGI and seed seq first so they're unlikely to be filtered out:
cat ../$GENE.fasta $GENE.JGI.query.sorted.fn.dd.fasta $GENE.jackhmmer.sorted.fasta $GENE.NCBI.query.fn.dd.fasta >> $GENE.all_ref.raw.fasta

# cluster at some threshold to reduce redundancy. We lose the aligment doing this. This will be used for the tree.
$USEARCH_PATH/usearch -cluster_fast $GENE.all_ref.raw.fasta -id 0.$USEARCH_THRESH -centroids $GENE.all_ref.id"$USEARCH_THRESH".fasta -uc $GENE.all_ref.id"$USEARCH_THRESH".uc

# align
mafft $GENE.all_ref.id"$USEARCH_THRESH".fasta > $GENE.all_ref.id"$USEARCH_THRESH".aln.fasta

# make the tree
FastTree -log $GENE.tree.log $GENE.all_ref.id"$USEARCH_THRESH".aln.fasta > $GENE.tree

# more refpkg stuff
make_seq_info_all.py $GENE.all_ref.id"$USEARCH_THRESH".aln.fasta
taxit update_taxids -d $TAX_DB -o seq_info_all.updated.csv seq_info_all.csv
cat seq_info_all.updated.csv | awk -F, '{print $2}' | sort | uniq | sed '/tax_id/d' > tax_ids_all.txt
taxit taxtable -d $TAX_DB  -t tax_ids_all.txt -o taxa.csv

# create the reference package:
taxit create -l $GENE -P $GENE.refpkg \
--taxonomy taxa.csv \
--aln-fasta $GENE.all_ref.id"$USEARCH_THRESH".aln.fasta \
--seq-info seq_info_all.updated.csv \
--tree-stats $GENE.tree.log \
--tree-file $GENE.tree \
--no-reroot

# now that we've recruited with marine ref, we'll build an aligment using our curated alignment.
hmmbuild $GENE.ref.hmm $GENE.all_ref.id"$USEARCH_THRESH".aln.fasta

cd ../hmmer_env

REFPKG="../hmmer_ref/$GENE.refpkg"
function pplacer_run {
	hmmalign -o $GENE.$1.aln.sto --mapali ../hmmer_ref/$GENE.all_ref.id"$USEARCH_THRESH".aln.fasta ../hmmer_ref/$GENE.ref.hmm $GENE."$1".query.sto
	seqmagick convert $GENE.$1.aln.sto $GENE.$1.aln.fasta
	$PPLACER_PATH/pplacer -c $REFPKG --keep-at-most 1 $GENE.$1.aln.fasta
}

for subject in $(cat $SUBJECT_LIST); do
  pplacer_run $subject
done

tar czf $GENE.gradients1.jplace.tar.gz *jplace
