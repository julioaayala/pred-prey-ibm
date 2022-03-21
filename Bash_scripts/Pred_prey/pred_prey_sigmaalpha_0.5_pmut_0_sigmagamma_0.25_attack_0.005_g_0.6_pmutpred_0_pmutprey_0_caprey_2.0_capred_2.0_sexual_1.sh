#!/bin/bash
#SBATCH -A lu2021-2-124
#SBATCH -N 4
#SBATCH -t 4:00:00
#SBATCH --mail-user=ju7141ay-s@student.lu.se
#SBATCH --mail-type=FAIL
#SBATCH -J pred_prey_0.5_0_0.25_0.005_0.6_0_3
#SBATCH -o Logs/pred_prey_0.5_0_0.25_0.005_0.6_0_3.out
#SBATCH -e Logs/pred_prey_0.5_0_0.25_0.005_0.6_0_3.err
#SBATCH -D /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm

module load matlab

cd /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm
matlab -r "main_pred_prey(0.5,0.25,0.005,0.6,3,0,0,1,2.0,2.0)"
