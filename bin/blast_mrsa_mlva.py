import argparse, os, glob, subprocess, textwrap
from pathlib import Path
from termcolor import colored

def getmylogo(pth):
    exec_globals = {}
    with open(pth, 'r') as lfile:
        exec(lfile.read(), exec_globals)
    logo = exec_globals.get('logo', None)
    return logo

def parse_arguments(logo):
    arg = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent(f"""
        {colored(logo, 'red', attrs=["bold"])}
        {colored('In silico MLVA typing for MRSA:', 'white', attrs=["bold", "underline"])}

        Blasts for all primers and repeat sequences.
        Only the repeat sequence of VNTR63_01 is used in downstream analysis.
-----------------------------------------------------------------------------------
        {colored('Example usage:', 'green', attrs=["bold", "underline"])}
            python {os.path.abspath(__file__)} 
            --input example/input_fasta/
            --output example/output_blastn/)
-----------------------------------------------------------------------------------
        """))
    arg.add_argument("-i", 
                    "--input", 
                    metavar="Path", 
                    help="Input directory with assembled fasta file", 
                    type=str, 
                    required=True)
    arg.add_argument("-o", 
                    "--output", 
                    metavar="Path", 
                    help="Output directory where you want to copy to", 
                    type=str, 
                    required=False)
    arg.add_argument("-pi", 
                    "--perc_identity",
                    metavar="INT", 
                    help="Percentage of identity to use to search for in the primers", 
                    type=parse_percentage, 
                    default=50, 
                    required=False)

    return arg.parse_args()

def determine_outdir(flg_out):
    if flg_out == None:
        outdir = os.path.abspath('') + '/output_blastn'
        Path(os.path.abspath(f"{outdir}")).mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {outdir}")
    else:
        Path(os.path.abspath(flg_out)).mkdir(parents=True, exist_ok=True)
        outdir = os.path.abspath(flg_out)
        print(f"Output directory: {outdir}")
    return outdir

def parse_percentage(value):
    value = int(value)
    if 0 <= value <= 100:
        return value
    else:
        raise argparse.ArgumentTypeError("Percentage must be in the range 0-100")

def main():
    current_file_path = os.path.abspath(__file__)
    parent_dir_path = os.path.dirname(os.path.dirname(current_file_path))
    primer_file = os.path.join(parent_dir_path, "files", "mrsa_mlva_primers.fasta")
    sequence_file = os.path.join(parent_dir_path, "files", "mrsa_mlva_sequenties.fasta")
    logo_path = os.path.join(parent_dir_path, "files", "logo.txt")
    flags = parse_arguments(getmylogo(logo_path))
    outdir = determine_outdir(flags.output)

    list_of_files = glob.glob(os.path.abspath(f"{flags.input}/*"))

    for key in list_of_files:
        basename = os.path.splitext(os.path.basename(key))[0]
        outputname = f"{outdir}/{basename}"
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
