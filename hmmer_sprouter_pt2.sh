# hmmer_sprouter_pt2.sh

GENE=$1
BASE_DIR=/mnt/nfs/ryan/Gradients1/runs/$GENE
cd $BASE_DIR

# Number of cores: 16 for gross and 8 for match
NCORES=$2
BITSCORE_CUTOFF=$3
USEARCH_THRESH=$4
PPLACER_PATH="/mnt/nfs/home/rgrous83/bin/pplacer"
SUBJECT_DIR="/mnt/nfs/ryan/Gradients1/mORFeus_v2"
SUBJECT_LIST="/mnt/nfs/ryan/Gradients1/mORFeus_v2/gradients1_morfeus_handles.txt"

# Re-situation directories
mkdir hmmer_env; cd hmmer_env

ref_hmm=$BASE_DIR/hmmer_ref/$GENE.ref.hmm
REFPKG="$BASE_DIR/$GENE.refpkg"

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

# now summarize the count information by treecolor groups:
# this will count each placement within a definined group (e.g., 'Bacillariophyta')
cd ../summed_csv
TREECOLORS="/mnt/nfs/home/rgrous83/scripts/treecolor/MarineRef2_plus_internal/treecolors_w_proks2.csv"
NORM_FACTORS="/mnt/nfs/ryan/Gradients1/gradients1.norm_factor_SUMS.csv"
count_pplacer_csv_by_taxonomy.py -eg -c $TREECOLORS -n $NORM_FACTORS $(ls $GENE.*.taxID.csv | head -1) > $GENE.g1.counts_results.csv # make a header
for csv in $(ls $GENE.*.taxID.csv); do count_pplacer_csv_by_taxonomy.py  -g -c $TREECOLORS $csv -n $NORM_FACTORS >> $GENE.g1.counts_results.csv; done

# normalize the raw counts by counts of 13 ribosomal proteins (RPx13):
RPx13="/mnt/nfs/ryan/Gradients1/RP/RPx13.g1.rp_avg.csv"
normalize_by_RP.py -n $BASE_DIR/summed_csv/$GENE.g1.counts_results.csv -d $RPx13 -o $BASE_DIR/$GENE.g1.counts_results.RPx13_norm.csv
