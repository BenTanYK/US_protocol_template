#!/bin/bash
#SBATCH -p main
#SBATCH -n 1
#SBATCH --gres=gpu:1 		
#SBATCH --mem=4G   
#SBATCH -o outfiles/%x-%j.out
#SBATCH -e outfiles/%x-%j.err

# Make directory to store slurm output
mkdir -p outfiles

# Initialize variables to store the numbers
RNUMBER=0

# Use flags r and j to store the specified values
while getopts "r:" opt; do
  case $opt in
    r)
      RNUMBER=$OPTARG
      ;;
    \?)
      echo "Usage: $0 -r <number> "
      exit 1
      ;;
  esac
done

# Shift off the processed options
shift $((OPTIND -1))

source /home/btan/miniconda3/etc/profile.d/conda.sh
conda activate openbiosim

python runMD.py $RNUMBER