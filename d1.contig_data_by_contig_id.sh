# d1.contig_data_by_contig_id.sh

# given a list of contig IDs, will extract information on PFAM annotations and counts:
# best to run from /mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99/diamond/combined/by_group/

# this is a list of contig_ids passed along on the command line:
group_name=$1
contig_ids=$1.contig_ids.txt
WORKING_DIR=/mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99

# Go into the group-specific directory:
cd $group_name

# first, pull out the matching contig_ids from the hmm_out.tab file for each collection of times:
for run in 0200 0600 1000 1400 1800 2200; do
echo "Retrieving PFAM annotations from $run"
grep -Ff $contig_ids $WORKING_DIR/pfam/diel1.bf100_id99."$run"_vs_PFAM.hmm_out.tab >> $group_name.d1.pfam.hmm_out.tab
done

# and also, lets get some normalized count information:
# convert contig_ids to nt ids (kallist counts are on the nt contigs):
cat $contig_ids | sed -r 's/_[0-9]+//g' > $group_name.nt_contig_ids.txt
# grab the header from one file (they are all the same)
run=0200; head -n 1 $WORKING_DIR/kallisto/"$run"_counts_full/combined/diel1.$run.bf100_id99.normed_counts.csv >> $group_name.d1.normed_counts.csv
for run in 0200 0600b 1000 1400 1800 2200; do
echo "Retrieving counts from $run"
grep -Ff $group_name.nt_contig_ids.txt $WORKING_DIR/kallisto/"$run"_counts_full/combined/diel1.$run.bf100_id99.normed_counts.csv >> $group_name.d1.normed_counts.csv
done

# and we'll also grab information on significant diel periodicity:
for peaktime in 4 8 12 16 20 24; do
echo "Retrieving peak times from $run"
grep -Ff $group_name.nt_contig_ids.txt /mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99/kallisto/RAIN_checked/$peaktime.fastaheaders.diel1.all.bf100_id99.filtered_counts.csv.rain.BH.sign.counts.txt > $group_name.d1.$peaktime.RAIN_sig.txt
done

cd ..
