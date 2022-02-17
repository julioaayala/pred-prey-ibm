#!/bin/bash
#SBATCH -A lu2021-2-124
#SBATCH -N 4
#SBATCH -t 4:00:00
#SBATCH --mail-user=ju7141ay-s@student.lu.se
#SBATCH --mail-type=FAIL
#SBATCH -J pred_prey_alpha_pmutprey_gamma_aP_efficiency_pmutpred_morphsinit
#SBATCH -o Logs/pred_prey_alpha_pmutprey_gamma_aP_efficiency_pmutpred_morphsinit.out
#SBATCH -e Logs/pred_prey_alpha_pmutprey_gamma_aP_efficiency_pmutpred_morphsinit.err
#SBATCH -D /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm

module load matlab

cd /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm
matlab -r "main_prey_pred(alpha,gamma,aP,efficiency,morphsinit,pmutprey,pmutpred)"