### In silico MLVA for MRSA script
This script will in-silico determine the MLVA for MRSA isolates based on a long-read only assembly.
It does so by blasting for all the primers and then determine if a product could have been formed that is <1200 bp.
For VNTR63_01 however the reverse primer was not found for 35% of tested isolates. So the approach here is to look for the repeat sequences that are within range of the forward primer
This makes it a different approach because the product size is not taken into account here as is the case with the in vitro approach.

### Example usage ###

bash run_pipeline.sh --input example/input_fasta/ --output example/output/