#!/bin/bash
#SBATCH -A lu2021-2-124
#SBATCH -N 4
#SBATCH -t 4:00:00
#SBATCH --mail-user=ju7141ay-s@student.lu.se
#SBATCH --mail-type=FAIL
#SBATCH -J pred_prey_0.5_0_0.75_0.015_0.8_0_3
#SBATCH -o Logs/pred_prey_0.5_0_0.75_0.015_0.8_0_3.out
#SBATCH -e Logs/pred_prey_0.5_0_0.75_0.015_0.8_0_3.err
#SBATCH -D /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm

module load matlab

cd /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm
matlab -r "main_pred_prey(0.5,0.75,0.015,0.8,3,0,0,0,1.0,1.0)"
