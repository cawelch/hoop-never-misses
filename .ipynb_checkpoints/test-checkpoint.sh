#!/bin/sh
echo "Starting SLURM job..."
/usr/bin/hostname

#SBATCH --job-name=test                     # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cawelch@davidson.edu    # Where to send mail
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=1gb                           # Job memory request
#SBATCH --time=06:00:00                     # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log         # Standard output and error log
pwd; hostname; date

#module load python

echo "Running plot script on a single CPU core"

python monte_carlo.py

date
