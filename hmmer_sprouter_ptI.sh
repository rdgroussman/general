### hmmer_sprouter_ptI.sh ###

GENE=$1
BASE_DIR=/mnt/nfs/ryan/Gradients1/runs/$GENE
mkdir $BASE_DIR; cd $BASE_DIR

REFSEQS_PATH="/mnt/nfs/ryan/MarineRef/with_virus/MarRef_w_virus.fn.fa"
# Number of cores: 16 for gross and 8 for match
NCORES=$2
BITSCORE_CUTOFF=$3
USEARCH_THRESH=$4
USEARCH_PATH="/mnt/nfs/home/rgrous83/bin/usearch"

# the seed hmm_profile should have this format: $GENE.hmm
hmm_profile=$BASE_DIR/$GENE.hmm

# Recruit reference sequences using hmm profile:
function hmmer_time_ref {
mkdir hmmer_ref; cd $BASE_DIR/hmmer_ref
GENE=$1
hmmsearch -T $BITSCORE_CUTOFF --incT $BITSCORE_CUTOFF --cpu $NCORES --tblout $GENE.MarRef2Plus.hmm_out.tab -A "$GENE".MarRef2Plus.query.sto $hmm_profile $REFSEQS_PATH
seqmagick convert $GENE.MarRef2Plus.query.sto $GENE.MarRef2Plus.query.fasta
# catch and remove bad terms in $GENE.MarRef2Plus.query.fasta:
sed -i "s/'//g" $GENE.MarRef2Plus.query.fasta
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

# convert the Newick tree to XML:
java -cp ~/bin/forester.jar org.forester.application.phyloxml_converter -f=nn -m hmmer_ref/$GENE.tree $GENE.ref_tree.xml
# color the tree by phylogeny:
TREECOLORS="treecolors_w_proks2.csv"
treecolor2.py -a -c $TREECOLORS -o $GENE.all.ref_tree.col.xml $GENE.ref_tree.xml
