
### hmmer_sprouter ###


GENE=$1
BASE_DIR=/mnt/nfs/ryan/diel1/runs/$GENE
mkdir $BASE_DIR; cd $BASE_DIR

REFSEQS_PATH="/mnt/nfs/ryan/MarineRef/with_virus/MarRef_w_virus.fn.fa"
# Number of cores: 16 for gross and 8 for match
NCORES=$2
BITSCORE_CUTOFF=$3
USEARCH_THRESH=$4
USEARCH_PATH="/mnt/nfs/home/rgrous83/bin/usearch"
PPLACER_PATH="/mnt/nfs/home/rgrous83/bin/pplacer"
SUBJECT_DIR="/mnt/nfs/ryan/diel1/mORFeus_v2"
SUBJECT_LIST="/mnt/nfs/ryan/diel1/mORFeus_v2/all_6tr_orfs40_handles.txt"

# the seed hmm_profile should have this format: $GENE.hmm
hmm_profile=$BASE_DIR/$GENE.hmm

# Recruit reference sequences using hmm profile:
function hmmer_time_ref {
mkdir hmmer_ref; cd $BASE_DIR/hmmer_ref
GENE=$1
hmmsearch -T $BITSCORE_CUTOFF --incT $BITSCORE_CUTOFF --cpu $NCORES --tblout $GENE.MarRef2Plus.hmm_out.tab -A "$GENE".MarRef2Plus.query.sto $hmm_profile $REFSEQS_PATH
seqmagick convert $GENE.MarRef2Plus.query.sto $GENE.MarRef2Plus.query.fasta
$USEARCH_PATH/usearch -sortbylength $GENE.MarRef2Plus.query.fasta -fastaout $GENE.MarRef2Plus.lensort.fasta
$USEARCH_PATH/usearch -cluster_fast $GENE.MarRef2Plus.lensort.fasta -id 0.$USEARCH_THRESH -centroids $GENE.MarRef2Plus.id$USEARCH_THRESH.fasta -uc $GENE.MarRef2Plus.id$USEARCH_THRESH.uc
hmmalign -o $GENE.MarRef2Plus.id$USEARCH_THRESH.aln.sto $hmm_profile $GENE.MarRef2Plus.id$USEARCH_THRESH.fasta
seqmagick convert $GENE.MarRef2Plus.id$USEARCH_THRESH.aln.sto $GENE.MarRef2Plus.id$USEARCH_THRESH.aln.fasta

hmmbuild $GENE.ref.hmm $GENE.MarRef2Plus.id$USEARCH_THRESH.aln.sto
ref_hmm=$BASE_DIR/hmmer_ref/$GENE.ref.hmm
# make the tree
FastTree -log $GENE.tree.log $GENE.MarRef2Plus.id$USEARCH_THRESH.aln.fasta > $GENE.tree
}

# collect and process the reference information:
hmmer_time_ref $GENE

# Process the reference sequence information for the refpkg:
TAX_DB="/mnt/nfs/home/rgrous83/NCBI/taxonomy.db"
make_seq_info_all.py $GENE.MarRef2Plus.id$USEARCH_THRESH.aln.fasta
taxit update_taxids -d $TAX_DB -o seq_info_all.updated.csv seq_info_all.csv
cat seq_info_all.updated.csv | awk -F, '{print $2}' | sort | uniq | sed '/tax_id/d' > tax_ids_all.txt
taxit taxtable -d $TAX_DB  -t tax_ids_all.txt -o taxa.csv

cd $BASE_DIR
# Create the refpkg:
taxit create -l $GENE -P $GENE.refpkg \
--taxonomy hmmer_ref/taxa.csv \
--aln-fasta hmmer_ref/$GENE.MarRef2Plus.id$USEARCH_THRESH.aln.fasta \
--seq-info hmmer_ref/seq_info_all.updated.csv \
--tree-stats hmmer_ref/$GENE.tree.log \
--tree-file hmmer_ref/$GENE.tree \
--no-reroot

# build the hmm profile:
REFPKG="$BASE_DIR/$GENE.refpkg"

# Re-situation directories
mkdir hmmer_env; cd hmmer_env

# Use the hmm profile to search against the environmental data; align to the profile and run pplacer
function hmmer_time_env {
SUBJECT_FASTA="$1".6tr.orfs40.fasta.gz
hmmsearch -T $BITSCORE_CUTOFF --incT $BITSCORE_CUTOFF --cpu $NCORES --tblout $GENE.$1.hmm_out.tab -A $GENE.$1.query.sto $ref_hmm $SUBJECT_DIR/$SUBJECT_FASTA
hmmalign -o $GENE.$1.aln.sto --mapali ../hmmer_ref/$GENE.MarRef2Plus.id$USEARCH_THRESH.aln.fasta $ref_hmm $GENE."$1".query.sto
seqmagick convert $GENE.$1.aln.sto $GENE.$1.aln.fasta
$PPLACER_PATH/pplacer -c $REFPKG --keep-at-most 1 $GENE.$1.aln.fasta
}

# run hmmer & pplacer against environmental seqs:
for subject in $(cat $SUBJECT_LIST); do
hmmer_time_env $subject
done

# getting a sample list (bloom)
cd $BASE_DIR/hmmer_env
ls $GENE.*jplace | awk -F. {'print $2'} | sort | uniq > ../sample_list.txt
mkdir ../summed_csv

# process CSVs by summed samples:
for sample in $(cat ../sample_list.txt); do
guppy to_csv -o ../summed_csv/$GENE.$sample.taxID.csv $GENE.$sample.*.jplace
done

cd ../summed_csv
TREECOLORS="/mnt/nfs/home/rgrous83/scripts/treecolor/MarineRef2_plus_internal/treecolors_w_proks.csv"
NORM_FACTORS="/mnt/nfs/ryan/diel1/diel1_NORM_FACTORS_v2.csv"
count_pplacer_csv_by_taxonomy.py -eg -c $TREECOLORS -n $NORM_FACTORS $(ls $GENE.*.taxID.csv | head -1) > $GENE.d1.counts_results.csv # make a header
for csv in $(ls $GENE.*.taxID.csv); do count_pplacer_csv_by_taxonomy.py  -g -c $TREECOLORS $csv -n $NORM_FACTORS >> $GENE.d1.counts_results.csv; done
