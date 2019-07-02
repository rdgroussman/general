# d1.contig_data_by_pfam.sh
# best to run from /mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99/extracted_by_pfam

# this is our pfam_id passed along on the command line:
pfam_id=$1
# make a subdirectory for the pfam id:
mkdir $pfam_id; cd $pfam_id
echo "Retrieving data for PFAM $pfam_id"
# all of the data sources are here:
WORKING_DIR=/mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99

# first, pull out the matching pfam_id from the hmm_out.tab file for each collection of times:
echo "Retrieving PFAM annotations form $pfam_id"
grep $pfam_id $WORKING_DIR/pfam/diel1.best_pfam_combined.deR.csv >> $pfam_id.d1.best_pf.hmm_out.csv

# get the taxids from this file:
cat $pfam_id.d1.best_pf.hmm_out.csv | awk -F"," {'print $1'} > $pfam_id.d1.contig_ids.txt

# now, get the corresponding tax mappings from the diamond runs:
for run in 0200 0600 1000 1400 1800 2200; do
echo "Retrieving DIAMOND annotations from $run"
grep -Ff $pfam_id.d1.contig_ids.txt $WORKING_DIR/diamond/d1.all_"$run".bf100.id99.vs_MarRef.lca.tab >> $pfam_id.d1.lca_tax.tab
done

# and also, lets get some normalized count information:
# convert contig_ids to nt ids (kallist counts are on the nt contigs):
cat $pfam_id.d1.contig_ids.txt | sed -r 's/_[0-9]+//g' > $pfam_id.d1.nt_contig_ids.txt
# grab the header from one file (they are all the same)
run=0200; head -n 1 $WORKING_DIR/kallisto/"$run"_counts_full/combined/diel1.$run.bf100_id99.normed_counts.csv >> $pfam_id.d1.normed_counts.csv
for run in 0200 0600b 1000 1400 1800 2200; do
echo "Retrieving normed_counts counts from $run"
grep -Ff $pfam_id.d1.nt_contig_ids.txt $WORKING_DIR/kallisto/"$run"_counts_full/combined/diel1.$run.bf100_id99.normed_counts.csv >> $pfam_id.d1.normed_counts.csv
done

# For the F in FPKMx:
# collect from /mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99/kallisto/*_counts_full/combined/diel1.*.kallisto_est_counts.tsv
# header:
run=0200; head -n 1 $WORKING_DIR/kallisto/"$run"_counts_full/combined/diel1.$run.bf100_id99.kallisto_est_counts.tsv >> $pfam_id.d1.kallisto_est_counts.tsv
for run in 0200 0600b 1000 1400 1800 2200; do
echo "Retrieving est_counts from $run"
grep -Ff $pfam_id.d1.nt_contig_ids.txt $WORKING_DIR/kallisto/"$run"_counts_full/combined/diel1.$run.bf100_id99.kallisto_est_counts.tsv >> $pfam_id.d1.kallisto_est_counts.tsv
done

# for the K in FPKMx:
# grab the lengths from the length file:
echo "target_id,length" > $pfam_id.d1.contig_lengths.csv
grep -Ff $pfam_id.d1.nt_contig_ids.txt $WORKING_DIR/kallisto/combined/diel1.bf100_id99.contig_lengths.csv >> $pfam_id.d1.contig_lengths.csv


# get RAIN data:
for run in 0200 0600b 1000 1400 1800 2200; do
echo "Retrieving RAIN results from $run"
grep -Ff $pfam_id.d1.nt_contig_ids.txt /mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99/kallisto/RAIN_checked/diel1.$run.bf100_id99.filtered_counts.csv.rain.BH.csv >> $pfam_id.d1.rain.csv
done

# and go back to the previous directory!
cd ..
