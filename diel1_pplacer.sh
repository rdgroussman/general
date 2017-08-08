



## runs a sample against the 192 samples ##
HANDLE=$1
SUBJECT_DIR="/mnt/ryan/diel1/6tr"
# SUBJECT_LIST="/mnt/ryan/diel1/6tr_handles.list" # the real one
# COLLAPSED_SUBJECT_LIST=""

### TESTING PURPOSES ###
SUBJECT_LIST="/mnt/ryan/diel1/test_6tr_handles.list" # TESTING
COLLAPSED_SUBJECT_LIST="test_6tr_collapsed_handles.list" #TESTING
########################


# Do an HMM search on synthetic Thaps transcriptome and get output in .sto format
for subject in $(cat $SUBJECT_LIST); do
  # build the filepath for subject files
  SUBJECT_FILE="$subject"_combined.6tr.fasta.gz

  # do the hmmsearch on each
  hmmsearch -A "$HANDLE"."$subject".query.sto "$HANDLE".hmm $SUBJECT_DIR/$SUBJECT_FILE
  # output: "$HANDLE"."$subject".query.sto
  # options for modification: change E-value threshold
  # -E <x>  where x > 0
  # use parallel processing: --cpu 16

  # use hmmalign to align query hits to the reference alignment
  hmmalign -o "$HANDLE"."$subject".aln.sto --mapali "$HANDLE".Aln.sto "$HANDLE".hmm "$HANDLE"."$subject".query.sto
  # output: "$HANDLE"."$subject".aln.sto

  # run pplacer using refpkg
  pplacer -c "$HANDLE".refpkg --keep-at-most 1 "$HANDLE"."$subject".aln.sto
  # output:

  # we are encountering an error in most of our searches here:
  #   Running pplacer v1.1.alpha14-1-gc47d954 analysis on LOV.root.S14C1_C_2200.H2NVM.aln.sto...
  # Found reference sequences in given alignment file. Using those for reference alignment.
  # Pre-masking sequences... Uncaught exception: Failure("* is not a known base in 561|2768_5/1-66")
  # Fatal error: exception Failure("* is not a known base in 561|2768_5/1-66")

  # so we should go with our mORFeus seqs then.

done

# did test run with HANDLE="LOV.root" in /mnt/ryan/diel1/runs/LOV/
# to do: look at time between each file to see how long it all takes.

# merge the machine-run results and export to csv
for subject in $(cat $COLLAPSED_SUBJECT_LIST); do
  guppy merge --split-csv -o "$subject"_combined.jplace "$subject"*.jplace
  guppy to_csv -o "$HANDLE"."$subject".TaxID.aln.csv "$HANDLE"."$subject".aln.jplace
done
