
## Build sh scripts for prey
for alpha in $(cat Params/sigma_alpha.txt);
  do for c_a_prey in $(cat Params/prey_choosiness_a.txt);
    do for sexual_rep in 0 1;
      do sed "s/alpha/$alpha/g" prey.sh | sed "s/c_a_prey/$c_a_prey/g" | sed "s/morphsinit/1/g" | sed "s/pmutprey/1e-5/g" | sed "s/is_sexual/$sexual_rep/g" \
      > Prey/prey_sigmaalpha_$(echo $alpha)_pmut_1e-5_ca_$(echo $c_a_prey)_sexual_$(echo $sexual_rep).sh;
    done;
  done;
done

## TODO: "Build"
## Build sh scripts for pred_prey (3 prey morphs, non-evolving)
for gamma in $(cat Params/sigma_gamma.txt);
  do for aP in $(cat Params/pred_attack.txt);
    do for efficiency in $(cat Params/pred_efficiency.txt);
      do for c_a_pred in $(cat Params/pred_choosiness_a.txt);
        do for sexual_rep in 0 1;
          do sed "s/alpha/0.5/g" pred_prey.sh| sed "s/gamma/$gamma/g" | sed "s/aP/$aP/g" | sed "s/efficiency/$efficiency/g" | \
           sed "s/morphsinit/3/g" | sed "s/morphsinit/3/g" | sed "s/pmutprey/0/g" | sed "s/pmutpred/0/g" | sed "s/is_sexual/$sexual_rep/g" |\
           sed "s/c_a_prey/$c_a_pred/g" | sed "s/c_a_pred/$c_a_pred/g" \
           > Pred_prey/pred_prey_sigmaalpha_0.5_pmut_0_sigmagamma_$(echo $gamma)_attack_$(echo $aP)_g_$(echo $efficiency)_pmutpred_0_pmutprey_0_caprey_$(echo $c_a_pred)_capred_$(echo $c_a_pred)_sexual_$(echo $sexual_rep).sh;
        done;
      done;
    done;
  done;
done

## Build sh scripts for pred_prey (3 prey morphs, evolving)
for gamma in $(cat Params/sigma_gamma.txt);
  do for aP in $(cat Params/pred_attack.txt);
    do for efficiency in $(cat Params/pred_efficiency.txt);
      do sed "s/gamma/$gamma/g" pred_prey.sh | sed "s/aP/$aP/g" | sed "s/efficiency/$efficiency/g" | \
      sed "s/alpha/0.5/g" | sed "s/morphsinit/3/g" | sed "s/morphsinit/3/g" | sed "s/pmutprey/0/g" | sed "s/pmutpred/1e-5/g" \
      > Pred_prey_evol/pred_prey_sigmaalpha_0.5_pmutprey_0_sigmagamma_$(echo $gamma)_attack_$(echo $aP)_g_$(echo $efficiency)_pmutpred_1e-5.sh;
    done;
  done;
done


## 1 morph
## Build sh scripts for pred_prey (1 prey morph, non-evolving)
for gamma in $(cat Params/sigma_gamma.txt);
  do for aP in $(cat Params/pred_attack.txt);
    do for efficiency in $(cat Params/pred_efficiency.txt);
      do sed "s/gamma/$gamma/g" pred_prey.sh | sed "s/aP/$aP/g" | sed "s/efficiency/$efficiency/g" | \
      sed "s/alpha/0.5/g" | sed "s/morphsinit/1/g" | sed "s/pmutprey/0/g" | sed "s/pmutpred/0/g" \
      > Pred_prey_1morph/pred_prey_sigmaalpha_0.5_pmut_0_sigmagamma_$(echo $gamma)_attack_$(echo $aP)_g_$(echo $efficiency)_pmutpred_0_pmutprey_0.sh;
    done;
  done;
done

## Build sh scripts for pred_prey (1 prey morph, evolving)
for gamma in $(cat Params/sigma_gamma.txt);
  do for aP in $(cat Params/pred_attack.txt);
    do for efficiency in $(cat Params/pred_efficiency.txt);
      do sed "s/gamma/$gamma/g" pred_prey.sh | sed "s/aP/$aP/g" | sed "s/efficiency/$efficiency/g" | \
      sed "s/alpha/0.5/g" | sed "s/morphsinit/1/g" | sed "s/pmutprey/1e-5/g" | sed "s/pmutpred/1e-5/g" \
      > Pred_prey_1morph_evol/pred_prey_sigmaalpha_0.5_pmutprey_1e-5_sigmagamma_$(echo $gamma)_attack_$(echo $aP)_g_$(echo $efficiency)_pmutpred_1e-5.sh;
    done;
  done;
done
