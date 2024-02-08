#!/bin/bash

set -x

##################
# in Silico MLVA #
##################

# EXAMPLE USAGE:
# bash run_pipeline.sh --input example/input_fasta --output example/output

version="In silico MLVA typer v0.1 (February 2024)"
dataecho=$(date +%y%m%d"_"%H"h"%M"m"%S"s")

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null 2>&1 && pwd )"
INPUT_CMD=""
OUTPUT_CMD=""
PATH_MASTER_YAML=$(echo "${DIR}/env/blastn_mlva.yaml")
MASTER_NAME=$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ')

function usage(){
	printf "Script usage :\n"
	printf "\t-h, --help				: Print this message and exit\n"
	printf "\t-v, --version				: Print the version and exit\n"
	printf "\t-i, --input				: Input directory with all your fasta files\n"
	printf "\t-o, --output			: Output directory, defaults to current dir + /output \n"
}

if [ $# == 0 ]
then
	usage
	exit 1
fi

while [[ $# -gt 0 ]]
do
    case "$1" in
    -v|--version)
        echo "$version"
        usage
        exit 0
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    -i|--input) 
        INPUT="$2";
        shift
        ;;
    -o|--output) 
        OUTPUT="$2";
        shift
        ;;    
    --) shift; break;;
    esac
    shift
done

# Determining the input directory
if [ -n "${INPUT}" ]
then
    export INPUT_DIR=$(realpath ${INPUT})
    INPUT_CMD="--input ${INPUT_DIR}"
else
    printf "No input was given\n"
    usage
    exit 1
fi

# Determining the output directory
if [ -n "${OUTPUT}" ]
then 
    export OUTPUT_DIR=$(realpath ${OUTPUT})
    OUTPUT_CMD=$(echo "--output ${OUTPUT_DIR}")
else
    printf "No output directory was given\n"
    usage
    exit 1
fi
mkdir -p "${OUTPUT_DIR}"/log/cluster


if [ ! -f $PATH_MASTER_YAML ]
    then
        printf "Warning! "$PATH_MASTER_YAML" file was not found! Can not run pipeline.\n"
        echo "Exiting. "
        exit 1
fi

set +ue # Turn bash strict mode off because that breaks conda

if [[ $PATH != *${MASTER_NAME}* ]]
    then # If echo $PATH does not contain the conda env name I will force install it.
    echo 'The environment '${MASTER_NAME}' is currently not in your PATH, will try to activate'; 
    source activate "${MASTER_NAME}" # Might be double..
    if ! source activate "${MASTER_NAME}"
    then # Only when it fails to activate this env it will force install.
        echo 'Could not find '${MASTER_NAME}' . Am now going to create a new environment with this name'
        # source activate mamba
        mamba env update -f ${PATH_MASTER_YAML}
        # conda env create -n "$MASTER_NAME" -f "$PATH_MASTER_YAML" # Old rule to install env but changed to mamba
        source activate "${MASTER_NAME}"
    fi
else # If it didn't fail it means I can activate it so am using that
    echo ${MASTER_NAME}' found, will now activate it.'; 
    source activate "${MASTER_NAME}"
fi

set -ue # Turn bash strict mode on again

python bin/blast_mrsa_mlva.py ${INPUT_CMD} ${OUTPUT_CMD}
sleep 5 # Wait for cluster to actually write the files
python bin/filter_mlva_blast.py ${INPUT_CMD} ${OUTPUT_CMD}
