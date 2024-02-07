import argparse
import os
import glob
import subprocess 
from pathlib import Path

# Example usage
# python /path/to/insilico_mlva/bin/blast_mrsa_mlva.py --input example/input_fasta/ --output example/output_blastn/ --perc_identity 50

def parse_percentage(value):
    value = int(value)
    if 0 <= value <= 100:
        return value
    else:
        raise argparse.ArgumentTypeError("Percentage must be in the range 0-100")

def main():
    arg = argparse.ArgumentParser(description="Blast MRSA MLVA")
    arg.add_argument("-i", 
                    "--input", 
                    metavar="Name", 
                    help="Input directory with assembled fasta file", 
                    type=str, 
                    required=True)
    arg.add_argument("-o", 
                    "--output", 
                    metavar="Name", 
                    help="Output directory where you want to copy to", 
                    type=str, 
                    required=False)
    arg.add_argument("-pi", 
                    "--perc_identity",
                    help="Percentage of identity to use to search for in the primers", 
                    type=parse_percentage, 
                    default=50, 
                    required=False)
    flags = arg.parse_args()

    if not flags.output:
        print("No output directory provided; using current directory + '/output' instead.")
        out = os.path.abspath('') + '/output_blastn'
    else:
        print("Output directory provided.")
        out = os.path.abspath(flags.output)

    Path(out).mkdir(parents=True, exist_ok=True)

    current_file_path = os.path.abspath(__file__)
    parent_dir_path = os.path.dirname(os.path.dirname(current_file_path))
    primer_file = os.path.join(parent_dir_path, "files", "mrsa_mlva_primers.fasta")
    sequence_file = os.path.join(parent_dir_path, "files", "mrsa_mlva_sequenties.fasta")

    list_of_files = glob.glob(os.path.abspath(f"{flags.input}/*"))

    for key in list_of_files:
        basename = os.path.splitext(os.path.basename(key))[0]
        outputname = f"{out}/{basename}"
        ### blast for primers:
        subprocess_cmd_primer_blast2input = f"bsub -q bio -n 1 -R 'rusage[mem=12G]' -R 'span[hosts=1]' -W 15 -M 16000 \
            \"blastn -query {key} \
            -subject {primer_file} \
            -out {outputname}_primers-blastn.csv \
            -word_size 7 \
            -perc_identity {flags.perc_identity} \
            -outfmt 10\""
        subprocess.Popen(subprocess_cmd_primer_blast2input, shell=True, stdout=subprocess.PIPE)
        ### blast for VNTR repeat sequences - used for VNTR63_01: 
        subprocess_cmd_repeat_blast2input = f"bsub -q bio -n 1 -R 'rusage[mem=12G]' -R 'span[hosts=1]' -W 15 -M 16000 \
            \"blastn -query {key} \
            -subject {sequence_file} \
            -out {outputname}_repeat-blastn.csv \
            -word_size 7 \
            -perc_identity {flags.perc_identity} \
            -outfmt 10\""
        subprocess.Popen(subprocess_cmd_repeat_blast2input, shell=True, stdout=subprocess.PIPE)

    print(f"Jobs sent to cluster for {len(list_of_files)} isolates.")

if __name__ == "__main__":
    main()
