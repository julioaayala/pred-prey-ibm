
## Build sh scripts for prey
for alpha in $(cat Params/sigma_alpha.txt); do sed "s/alpha/$alpha/g" prey.sh | sed "s/pmut/1e-4/g" > prey_sigmaalpha_$(echo $alpha)_pmut_1e-4.sh; done
for i in $(ls prey_*); do sbatch $i; done


## Build sh scripts for pred_prey (3 prey morphs, non-evolving)
for gamma in $(cat Params/sigma_gamma.txt);
  do for aP in $(cat Params/pred_attack.txt);
    do for efficiency in $(cat Params/pred_efficiency.txt);
      do sed "s/gamma/$gamma/g" pred_prey.sh | sed "s/aP/$aP/g" | sed "s/efficiency/$efficiency/g" | \
      sed "s/alpha/0.5/g" | sed "s/morphsinit/3/g" | sed "s/morphsinit/3/g" | sed "s/pmutprey/0/g" | sed "s/pmutpred/0/g" \
      > Pred_prey/pred_prey_sigmaalpha_0.5_pmut_0_sigmagamma_$(echo $gamma)_attack_$(echo $aP)_g_$(echo $efficiency)_pmutpred_0_pmutprey_0.sh;
    done;
  done;
done

## Build sh scripts for pred_prey (3 prey morphs, evolving)
for gamma in $(cat Params/sigma_gamma.txt);
  do for aP in $(cat Params/pred_attack.txt);
    do for efficiency in $(cat Params/pred_efficiency.txt);
      do sed "s/gamma/$gamma/g" pred_prey.sh | sed "s/aP/$aP/g" | sed "s/efficiency/$efficiency/g" | \
      sed "s/alpha/0.5/g" | sed "s/morphsinit/3/g" | sed "s/morphsinit/3/g" | sed "s/pmutprey/0/g" | sed "s/pmutpred/1e-4/g" \
      > Pred_prey/pred_prey_sigmaalpha_0.5_pmut_0_sigmagamma_$(echo $gamma)_attack_$(echo $aP)_g_$(echo $efficiency)_pmutpred_1e-4_pmutprey_0.sh;
    done;
  done;
done


# for i in $(ls prey_*); do sbatch $i; done
