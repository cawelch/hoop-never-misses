#!/bin/sh
echo "Starting SLURM job..."
/usr/bin/hostname

#SBATCH --job-name=monte_carlo_three        # Job name
#SBATCH --mail-user=cawelch@davidson.edu    # Where to send mail
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=2gb                           # Job memory request
#SBATCH --time=08:00:00                     # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log         # Standard output and error log
pwd; hostname; date

#module load python

echo "Running plot script on a single CPU core"

python monte_carlo.py > monte_carlo_output.txt 
sendmail cawelch@davidson.edu  < monte_carlo_output.txt


date
