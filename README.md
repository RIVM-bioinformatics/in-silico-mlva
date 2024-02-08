### In silico MLVA for MRSA script
This script will in-silico determine the MLVA for MRSA isolates based on a long-read only assembly.
It does so by blasting for all the primers and then determine if a product could have been formed that is <1200 bp.
For VNTR63_01 however the reverse primer was not found for 35% of tested isolates. So the approach here is to look for the repeat sequences that are within range of the forward primer
This makes it a different approach because the product size is not taken into account here as should be the case. 


#TODO It still crashes if it cant find the forward VNTR63 primer (local variable 'forward_pos1' referenced before assignment), in get_number_repeats

#TODO Check if the sizes found are at locations within range of their repeat sequences. Currently only done for VNTR63_01

#TODO if there are overlapping repeat sequences it cant tell the the right order, so it should only do this per contig found, then pick the one with the highest bitscore

### Example usage ###

bash run_pipeline.sh --input example/input_fasta/ --output example/output/