
### This script will in-silico determine the MLVA for MRSA isolates based on a long-read only assembly.
### It does so by blasting for all the primers and then determine if a product could have been formed that is <1200 bp.
### For VNTR63_01 however the reverse primer was not found for 35% of tested isolates. So the approach here is to look for the repeat sequences that are within range of the forward primer
### This makes it a different approach because the product size is not taken into account here as should be the case. 
###  to add both flanking regions + repeat sizes + primers sizes. But how to determine the flank region near the reverse primer? Can check in BN or other sequences if it's conserved.

#TODO It still crashes if it cant find the forward VNTR63 primer (local variable 'forward_pos1' referenced before assignment), in get_number_repeats

#TODO Check if the sizes found are at locations within range of their repeat sequences. Currently only done for VNTR63_01

#TODO if there are overlapping repeat sequences it cant tell the the right order, so it should only do this per contig found, then pick the one with the highest bitscore

inputdir="/path/to/insilico_mlva/example/input_fasta"
outputdir="/path/to/insilico_mlva/example"

### Example usage ###
### Activate your blastn environment and tools required for filter script ###

conda activate blast4mash
python /mnt/scratch_dir/landmanf/gitlabrivm/insilico_mlva/bin/blast_mrsa_mlva.py -i ${inputdir} -o ${outputdir}/output_blastn/

### wait for jobs ###

for file in ${inputdir}/*
do
base=$(basename "$file")
base_no_ext="${base%.fasta}"
primerfile="${outputdir}/${base_no_ext}_primers-blastn.csv"
repeatfile="${outputdir}/${base_no_ext}_repeat-blastn.csv"
python /path/to/insilico_mlva/bin/filter_mlva_blast.py -bp ${primerfile} -br ${repeatfile} -o  ${outputdir}/output_mlva_typing/
done
