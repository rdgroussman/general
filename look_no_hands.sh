

# Ryan Groussman
# Armbrust Lab 2017

# look_no_hands.sh


# Purpose:
# Will perform a 'quick and dirty' analysis of a given gene with the
# intention of providing quick visual results.
# If results look interesting, follow up with a careful, curated analysis.

# Declare the gene handle used during the analysis:
GENE=$1

# Collect the initial working directory
BASE_DIR=$PWD

# Number of cores. Match = 8, Gross = 16.
NCORES=16

# "97" will search at -id 0.97 and save $GENE.id97.fasta
USEARCH_THRESH=97

# Constant paths for databases, scripts, binaries, etc

# FASTA reference paths (local macbook):
# MARINEREF_PATH=""
# MMETSP_PATH="/Users/rgroussman/data/MMETSP/MMETSP.nr_clust.pep.fa"
# JGI_PATH="/Users/rgroussman/data/JGI/JGI_microeuks.aa.fasta"
# NCBI_PATH="/Users/rgroussman/data/NCBI/opisthokont_transcriptomes/marine_opisthokonts.fasta"

# # FASTA reference paths (bloom:match):
MMETSP_PATH="/mnt/nfs/ryan/MMETSP/MMETSP.nr_clust.pep.fa"
JGI_PATH="/mnt/nfs/ryan/JGI/JGI_microeuks.aa.fasta"
NCBI_PATH="/mnt/nfs/ryan/NCBI/extra_opisthokonts/marine_opisthokonts.fasta"
MARINEREF_PATH="/mnt/nfs/ryan/MarineRef/MarineRefII_ALL_SEQS.fa"
MARINEREF_INFO="/mnt/nfs/ryan/MarineRef/MarineRefII_seqIDinfo.csv"

# Paths for specific programs and scripts:
USEARCH_PATH="/mnt/nfs/home/rgrous83/bin/usearch"
PPLACER_PATH="/mnt/nfs/home/rgrous83/bin/pplacer"

# diel1 directory:
SUBJECT_DIR="/mnt/nfs/ryan/diel1/mORFeus/"
# SUBJECT_LIST="/mnt/nfs/ryan/diel1/mORFeus/all_6tr_orfs40_handles.txt" # full list of handles
# reduced version of full subject list from above:
# for run in "H5C5H" "H2NVM"; do grep "$run" $SUBJECT_LIST >> test.handles.txt; done
SUBJECT_LIST="/mnt/nfs/ryan/diel1/mORFeus/test.handles.txt"

TAX_DB="/mnt/nfs/home/rgrous83/NCBI/taxonomy.db"


# folder for reference hmmer files and environmental hmmer files:
mkdir hmmer_ref
mkdir hmmer_env

# Create an hmm profile for the initial fasta file:
# hmmbuild hmmer_ref/$GENE.hmm $GENE.aln.fasta

function hmmer_time {
	echo "Looking for $GENE in $2..."
	hmmsearch --cpu $NCORES --tblout hmmer_ref/$GENE.$1.hmm_out.tab --cut_tc -A hmmer_ref/$GENE.$1.ref.sto $GENE.hmm $2
}

hmmer_time JGI $JGI_PATH # very fast
hmmer_time NCBI_opistho $NCBI_PATH # very fast
hmmer_time MarRef $MARINEREF_PATH # takes a while...

# hmmer_time MMETSP $MMETSP_PATH

cd hmmer_ref
# convert to FASTA and sort by length; these steps are fast.
for sto in $(ls $GENE*ref.sto); do
	seqmagick convert $sto $sto.fasta
	$USEARCH_PATH/usearch -sortbylength $sto.fasta -fastaout $sto.sorted.fasta
done

# Rename the MarRefII seqs; takes a few moments.
python /mnt/nfs/home/rgrous83/scripts/Rename_MarineRefII_Seqs_fixed.py -I $MARINEREF_INFO -S $GENE.MarRef.ref.sto.sorted.fasta
fnames2.py -tmc $GENE.JGI.ref.sto.sorted.fasta
fnames2.py -tmc $GENE.NCBI_opistho.ref.sto.fasta # note that usearch drops the description so we lose tax_ids; so we're not usig the sorted ones.

# Add up all the reference sequences with JGI first so they're unlikely to be filtered out:
cat $GENE.JGI.ref.sto.sorted.fn.fasta named_"$GENE".MarRef.ref.sto.sorted.fasta $GENE.NCBI_opistho.ref.sto.fn.fasta >> $GENE.all_ref.raw.fasta

# cluster at some threshold to reduce redundancy:
$USEARCH_PATH/usearch -cluster_fast $GENE.all_ref.raw.fasta -id 0.$USEARCH_THRESH -centroids $GENE.ref.id"$USEARCH_THRESH".fasta -uc $GENE.ref.id"$USEARCH_THRESH".uc

# de-dupe
seqmagick mogrify --deduplicate-sequences $GENE.ref.id"$USEARCH_THRESH".fasta

# make the alignment, tree, and reference package:
# mafft $GENE.ref.id"$USEARCH_THRESH".fasta > $GENE.ref.id"$USEARCH_THRESH".aln.fasta

# OR we align with hmmalign - perhaps better.
hmmalign --trim -o $GENE.ref.aln.sto ../$GENE.hmm $GENE.ref.id"$USEARCH_THRESH".fasta
seqmagick convert $GENE.ref.aln.sto $GENE.ref.aln.fasta

make_seq_info_all.py $GENE.ref.id"$USEARCH_THRESH".aln.fasta

#####

# update taxids in seq_info_all.csv to seq_info_all.updated.csv
taxit update_taxids -d $TAX_DB -o seq_info_all.updated.csv seq_info_all.csv

### may need to do something there ###

cat seq_info_all.updated.csv | awk -F, '{print $2}' | sort | uniq | sed '/tax_id/d' > tax_ids_all.txt

taxit taxtable -d $TAX_DB  -t tax_ids_all.txt -o taxa.csv

FastTree -log $GENE.tree.log $GENE.ref.aln.fasta > $GENE.tree

# create the reference package:
taxit create -l $GENE -P $GENE.refpkg \
    --taxonomy taxa.csv \
    --aln-fasta $GENE.ref.aln.fasta \
    --seq-info seq_info_all.updated.csv \
    --tree-stats $GENE.tree.log \
    --tree-file $GENE.tree \
		--no-reroot

# now we'll pull out the seqs from the query set:
cd ../hmmer_env

function pplacer_run {
  # build the filepath for subject files
  # SUBJECT_FASTA="$1"_combined.6tr.orfs40.fasta.gz
  # hmmsearch --cpu $NCORES --tblout $GENE.$1.hmm_out.tab --cut_tc -A $GENE.$1.query.sto ../$GENE.hmm $SUBJECT_DIR/$SUBJECT_FASTA
	# Here's the step where we might insert filtering of STO output based on valuesin the .tab file.
	# convert output to FASTA and concatenate with the reference FASTA:
	# seqmagick convert $GENE.$1.query.sto $GENE.$1.query.fasta
	cat $GENE.$1.query.fasta ../hmmer_ref/$GENE.ref.aln.fasta >> $GENE.$1.ref_added.fasta
  # use hmmalign to align query hits to the reference alignment
  hmmalign --trim --informat FASTA -o $GENE.$1.aln.sto ../"$GENE".hmm $GENE.$1.ref_added.fasta
  # run pplacer using refpkg. might need to convert to FASTA first if it throws out problems:
	seqmagick convert $GENE.$1.aln.sto $GENE.$1.aln.fasta
  $PPLACER_PATH/pplacer -c ../hmmer_ref/$GENE.refpkg --keep-at-most 1 $GENE.$1.aln.fasta
}


for subject in $(cat $SUBJECT_LIST); do
  pplacer_run $subject
done
