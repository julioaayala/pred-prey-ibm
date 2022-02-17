#!/bin/bash
#SBATCH -A lu2021-2-124
#SBATCH -N 4
#SBATCH -t 4:00:00
#SBATCH --mail-user=ju7141ay-s@student.lu.se
#SBATCH --mail-type=FAIL
#SBATCH -J Prey_0.5_1e-4
#SBATCH -o Prey_0.5_1e-4.out
#SBATCH -e Prey_0.5_1e-4.err
#SBATCH -D /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm

module load matlab

cd /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm
matlab -r "main_prey(0.5,1e-4)"
