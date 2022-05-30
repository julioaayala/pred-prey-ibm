#!/bin/bash
#SBATCH -A lu2021-2-124
#SBATCH -N 4
#SBATCH -t 8:00:00
#SBATCH --mail-user=ju7141ay-s@student.lu.se
#SBATCH --mail-type=FAIL
#SBATCH -J Prey_alpha_pmutprey_morphsinit_is_sexual_c_a_prey.out
#SBATCH -o Prey_alpha_pmutprey_morphsinit_is_sexual_c_a_prey.out
#SBATCH -e Prey_alpha_pmutprey_morphsinit_is_sexual_c_a_prey.out
#SBATCH -D /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm

## --------------------------
## Julio Ayala
## February 2022
## prey.sh
## Template for execution of prey dynamics in Matlab, with a template from SNIC-LUNARC
## --------------------------

module load matlab

cd /home/ayalaj/Pred_prey_dynamics/pred-prey-ibm
matlab -r "main_pred_prey(alpha,morphsinit,pmutprey,is_sexual,c_a_prey)"
