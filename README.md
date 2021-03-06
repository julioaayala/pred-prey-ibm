# Pred-prey-ibm: A tool to model predator-prey eco-evolutionary feedbacks through an individual trait-based model
### Julio Ayala
#### May 2022
A repository with code associated project titled **"Diversification of predators in multi-trophic communities: A trait-based theoretical approach"**. This code is part of a master's thesis (60 cr) in bioinformatics in Lund University.

## Requirements
For the main tool for simulations:
- Matlab (Tested on vR2021b)
For analysis and plotting
- R (Tested on v4.0.3)
  - dplyr (v1.0.4)
  - ggplot2 (v3.3.3)
  - gridExtra (v2.3)

## Introduction
This repository contains 3 components:
- The main tool in Matlab, which contains the main code for simulations and is located in the root folder of this repository.
- Command line scripts located in the **Bash_scripts** folder, and contains templates for prey-only and predator-prey simulations. It also contains a folder named **Params**, which contain all the parameters the project was run on.
- A script in R, used for cluster analysis and plotting

## Execution
In order to avoid errors when executing the different scripts, file structure should be as follows:
```
.
├── Bash_scripts
│   ├── Params
│   │   ├── pmut_pred.txt
│   │   ├── pred_attack.txt
│   │   ├── pred_choosiness_a.txt
│   │   ├── pred_efficiency.txt
│   │   ├── prey_choosiness_a.txt
│   │   ├── sigma_alpha.txt
│   │   └── sigma_gamma.txt
│   ├── Pred_prey
│   ├── Pred_prey_1morph_evol
│   ├── Pred_prey_1morph_evol_sexual
│   ├── Pred_prey_evol_asexual
│   ├── Pred_prey_evol_sexual
│   ├── Prey
│   ├── builder.sh
│   ├── pred_prey.sh
│   └── prey.sh
├── Genetics.m
├── Individual.m
├── Plots
│   └── pred_prey_analysis.R
├── Population.m
├── Predator.m
├── Prey.m
├── Resource.m
├── Results
│   ├── Pred_prey_1morph
│   ├── Pred_prey_1morph_sexual
│   └── Prey
├── deterministic_main.m
├── fitness_landscape_preyonly.m
├── fitness_landscape_test.m
├── main_ode.m
├── main_pred_prey.m
└── reproduce.m
```

Empty folders for results and scripts can be created as:
``` bash
mkdir -p Results/{Pred_prey_1morph,Pred_prey_1morph_sexual,Prey} \
Bash_scripts/{Prey,Pred_prey,Pred_prey_evol_asexual,Pred_prey_evol_sexual,\
Pred_prey_1morph_evol,Pred_prey_1morph_evol_sexual}
````


### Main simulation
The main program can be executed from either the command line or Matlab by running:
```bash
matlab -r "main_pred_prey(alpha,gamma,aP,efficiency,morphsinit,pmutprey,pmutpred,is_sexual,c_a_prey,c_a_pred)"
```
Or, for prey-only dynamics:
```bash
matlab -r "main_pred_prey(alpha,morphsinit,pmutprey,is_sexual,c_a_prey)"
```

Where:
```
alpha = alpha niche width for prey
gamma = gamma niche width for predators
aP = Predator base attack rate
efficiency = Predator feeding efficiency
morphsinit = Initial number of prey morphs
pmutprey = Prey mutation rate
pmutpred = Predator mutation rate
is_sexual = Mode of reproduction: 1 for sexual, 0 for asexual
c_a_prey = Prey choosiness parameter
c_a_pred = Predator choosiness parameter
```


### Execution from bash script
Alternatively, scripts can be generated by running different snippets in `Bash_scripts/builder.sh`, which use parameter files located in `Bash_scripts/Params`. Parameter files contain line separated values:
```
pmut_pred.txt --> Predator mutation rates
pred_attack.txt --> Predator attack rates
pred_choosiness_a.txt --> Predator choosiness values
pred_efficiency.txt --> Predator feeding efficiency
prey_choosiness_a.txt --> Prey choosiness values
sigma_alpha.txt --> Prey resource niche width
sigma_gamma.txt --> Predator resource niche width
```

Execution of all scripts stores output in the `Results/` folder. Afterwards, they need to be moved manually to each of the different subfolders.

## Analysis and plotting
Analysis and plotting can be done by running the different snippets in `Plots/pred_prey_analysis.R`. First, the workspace needs to be set manually to the local `Results/` folder.
Plots are stored in .pdf files in the `Results/` folder.

---
For questions, bugs, or clarifications, contact ju7141ay-s[at]student.lu.se or ayala.j[at]live.com
