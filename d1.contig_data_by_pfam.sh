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
for run in 0200 0600 1000 1400 1800 2200; do
echo "Retrieving PFAM annotations from $run"
grep $pfam_id $WORKING_DIR/pfam/diel1.bf100_id99."$run"_vs_PFAM.hmm_out.tab >> $pfam_id.d1.hmm_out.tab
done

# get the taxids from this file:
cat $pfam_id.d1.hmm_out.tab | awk -F" " {'print $1'} > $pfam_id.d1.contig_ids.txt

# now, get the corresponding tax mappings from the diamond runs:
for run in 0200 0600 1000 1400 1800 2200; do
echo "Retrieving DIAMOND annotations from $run"
grep -Ff $pfam_id.d1.contig_ids.txt $WORKING_DIR/diamond/d1.all_"$run".bf100.id99.vs_MarRef.lca.tab >> $pfam_id.d1.lca_tax.tab
done

# and also, lets get some normalized count information:
# convert contig_ids to nt ids (kallist counts are on the nt contigs):
cat $pfam_id.d1.contig_ids.txt | sed -r 's/_[0-9]+//g' > $pfam_id.d1.nt_contig_ids.txt
for run in 0200 0600b 1000 1400 1800 2200; do
echo "Retrieving counts from $run"
grep -Ff $pfam_id.d1.nt_contig_ids.txt $WORKING_DIR/kallisto/"$run"_counts_full/combined/diel1.$run.bf100_id99.normed_counts.csv >> $pfam_id.d1.normed_counts.csv
done

# ok great! now let's turn to a python script to put it all together in one big file!
TAXA_CSV="/mnt/nfs/ryan/MarineRef/with_virus/taxa.csv"
combine_d1_contig_data.py $pfam_id -t $TAXA_CSV -o $pfam_id.d1.all_data.csv

# and go back to the previous directory!
cd ..
