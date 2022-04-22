#!/bin/bash
#SBATCH -A lu2021-2-124
#SBATCH -N 1
#SBATCH --tasks-per-node=8
#SBATCH -t 10:00:00
#SBATCH --mail-user=ju7141ay-s@student.lu.se
#SBATCH --mail-type=FAIL
#SBATCH -J pred_prey_alpha_pmutprey_gamma_aP_efficiency_pmutpred_morphsinit_c_a_prey_c_a_pred_is_sexual
#SBATCH -o Logs/pred_prey_alpha_pmutprey_gamma_aP_efficiency_pmutpred_morphsinit_c_a_prey_c_a_pred_is_sexual.out
#SBATCH -e Logs/pred_prey_alpha_pmutprey_gamma_aP_efficiency_pmutpred_morphsinit_c_a_prey_c_a_pred_is_sexual.err
#SBATCH -D /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm

module load matlab

cd /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm
matlab -r "main_pred_prey(alpha,gamma,aP,efficiency,morphsinit,pmutprey,pmutpred,is_sexual,c_a_prey,c_a_pred)"
