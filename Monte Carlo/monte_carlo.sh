#!/bin/bash
#SBATCH --job-name=monte_carlo              # Job name
#SBATCH --mail-user=cawelch@davidson.edu    # Where to send mail
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=10gb                          # Job memory request
#SBATCH --time=40:00:00                      # Time limit hrs:min:sec
#SBATCH --output=monte_carlo_%j.log         # Standard output and error log

echo "Starting SLURM job..."
/usr/bin/hostname

pwd; hostname; date

#module load python

#for i in {1...10}
#do
echo "Running plot script on a single CPU core"
python monte_carlo_two_d.py > monte_carlo.txt 
sendmail cawelch@davidson.edu  < monte_carlo.txt
#done
#wait
    


date
