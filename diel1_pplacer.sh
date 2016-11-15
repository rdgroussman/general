



## runs a sample against the 192 samples ##
HANDLE=$1
SUBJECT_DIR="/mnt/ryan/diel1/6tr"
# SUBJECT_LIST="/mnt/ryan/diel1/6tr_handles.list" # the real one
SUBJECT_LIST="/mnt/ryan/diel1/test_6tr_handles.list" # TESTING

# Do an HMM search on synthetic Thaps transcriptome and get output in .sto format
for subject in $(cat $SUBJECT_LIST); do
  # build the filepath for subject files
  SUBJECT_FILE="$subject"_combined.6tr.fasta.gz
  # do the hmmsearch on each
  hmmsearch -A "$HANDLE"."$subject".query.sto "$HANDLE".hmm $SUBJECT_DIR/$SUBJECT_FILE
  # use hmmalign to align query hits to the reference alignment
  hmmalign -o "$HANDLE"."$subject".aln.sto --mapali "$HANDLE".Aln.sto "$HANDLE".hmm "$HANDLE"."$subject".query.sto
  # run pplacer using refpkg
  pplacer -c "$HANDLE".refpkg --keep-at-most 1 "$HANDLE"."$subject".aln.sto

  # merge placefiles from the different machine runs:
  guppy merge -o "outfile.jplace" "placefiles"

  # export to .csv
  guppy to_csv -o "$HANDLE"."$subject".TaxID.aln.csv "$HANDLE"."$subject".aln.jplace
done
