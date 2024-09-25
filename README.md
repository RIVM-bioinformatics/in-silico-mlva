## Used for Genome Medicine paper
This pipeline was used to determine the in-silico MLVA for MRSA isolates for our 'Genomic surveillance of multidrug-resistant organisms based on long-read sequencing' submitted to Genome Medicine.

## In silico MLVA for MRSA script
This script will in-silico determine the MLVA for MRSA isolates based on a long-read only assembly.
It does so by blasting for all the primers and then determine if a product could have been formed that is <1200 bp.
For VNTR63_01 however the reverse primer was not found for 35% of tested isolates. So the approach here is to look for the repeat sequences that are within range of the forward primer
This makes it a different approach because the product size is not taken into account here as is the case with the in vitro approach.

## Usage 

bash run_pipeline.sh --input example/input_fasta/ --output example/output/

## Authors and acknowledgment
Pipeline written by Fabian Landman.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Project status
This pipeline has been uploaded for submission of our article only and shall not be further developed.