%------------------------------------------------------------
% Julio Ayala
% ju7141ay-s@student.lu.se
% May 2022
% Description: Main script to simulate predator-prey eco-evolutionary 
% dynamics. The script outputs a csv file in a ./Results folder located in
% the same location as the script.
% Usage:
% Run from terminal as: 
% For prey only:
%     matlab -r "main_pred_prey(alpha,morphsinit,pmutprey,is_sexual,c_a_prey)"
% For predator and prey:
%     matlab -r "main_pred_prey(alpha,gamma,aP,efficiency,morphsinit,pmutprey,pmutpred,is_sexual,c_a_prey,c_a_pred)"
% Where:
%     alpha = alpha niche width for prey
%     gamma = gamma niche width for predators
%     aP = Predator base attack rate
%     efficiency = Predator feeding efficiency
%     morphsinit = Initial number of prey morphs
%     pmutprey = Prey mutation rate
%     pmutpred = Predator mutation rate
%     is_sexual = Mode of reproduction: 1 for sexual, 0 for asexual
%     c_a_prey = Prey choosiness parameter
%     c_a_pred = Predator choosiness parameter
%------------------------------------------------------------

function main_pred_prey(varargin)
    %% Parameters ----------------------------------------------------
    % Constants
    t_end = 10000;              % Number of generations
    F = 2;                      % Fecundity rate
    K = 400;                    % Single resource system size (Related to carrying capacity)
    p_dispersal = 0;            % Probability of dispersal (Untested)
    N0N = int16(K);             % Initial prey pop size
    N0P = int16(K/20);          % Initial pred pop size
    bmax = 1;                   % Max attack rate
    a_0N = 2;                   % Prey base attack rate

    % Overriden parameters from terminal. Defaults when no parameters are
    % given
    p_mut_prey = 0;             % Probability of prey mutation
    p_mut_pred= 0;              % Probability of pred mutation
    sigma_alpha = 0.35;         % Prey-resource niche width
    sigma_gamma = 0.5;          % Predator niche width
    a_0P = 0.05;                % Predator attack rate
    g = 0.5;                    % Predator feeding efficiency
    num_populations = 1;        % One initial population, 1 for prey only; 2 for pred-prey
    morphs = 1;                 % Initial number of prey morphs
    sigma_beta = 0.5;           % Habitat niche width (Unused)
    num_resources = 3;          % Number of resources (Unused)
    num_habitats = 1;           % Habitats (Unused)
    c_a_pred = -log(0.5);       % Predator trait-choosiness
    c_a_prey = -log(0.5);       % Prey trait-choosiness
    c_ss_pred = 0;              % Predator sexual trait choosiness (Unused)
    c_ss_prey = 0;              % Prey sexual trait choosiness (Unused)
    is_sexual = 0;              % Mode of reproduction: 1 for sexual, 0 for asexual
    

    % Number of loci per trait
    % Prey
    loci_prey.alpha = 16;       % ecological alpha trait
    loci_prey.beta = 8;         % Niche beta trait (unused)
    loci_prey.dis = 8;          % Display trait (unused)
    loci_prey.pref = 8;         % Preference trait (unused)
    
    % Predator
    loci_pred.alpha = 32;       % ecological alpha(gamma) trait
    loci_pred.beta = 8;         % Niche beta trait (unused)
    loci_pred.dis = 8;          % Display trait (unused)
    loci_pred.pref = 8;         % Preference trait (unused)

    % Parameters from terminal
    if ~isempty(varargin)
        if length(varargin)==5  % Prey only
            num_populations = 1;
            sigma_alpha = varargin{1};
            morphs = varargin{2};
            p_mut_prey = varargin{3};
            is_sexual = varargin{4}; 
            c_a_prey =  varargin{5}; 
        else                    % Predator-prey
            num_populations = 2;
            sigma_alpha = varargin{1};
            sigma_gamma = varargin{2};
            a_0P = varargin{3};
            g = varargin{4};
            morphs = varargin{5};
            p_mut_prey = varargin{6};
            p_mut_pred = varargin{7};
            is_sexual = varargin{8};
            c_a_prey =  varargin{9};
            c_a_pred =  varargin{10};
        end
    end
    
    %% Initialization ----------------------------------------------------
    % Resources
    resource_abundance = zeros(t_end,num_habitats, num_resources);
    resources = Resource.empty;
    for i=1:num_habitats
        for j=1:num_resources
            resources(i,j) = Resource(j,i,K);
        end
    end

    % Individuals
    trait_frequency = {};   
    population_size = zeros(t_end,num_populations);
    population = Population.empty(num_populations,0);

    for i=1:num_populations
        if mod(i,2)==1                  % Prey
            if morphs>1                 % Initialize prey populations with more than 1 morph
                population(i) = Population("prey", 1:morphs, 1, a_0N, 1, sigma_alpha, sigma_beta, c_a_prey, c_ss_prey, N0N, p_mut_prey, loci_prey);
            else
                population(i) = Population("prey", 2, 1, a_0N, 1, sigma_alpha, sigma_beta, c_a_prey, c_ss_prey, N0N, p_mut_prey, loci_prey);
            end
        else                            % Predators
            population(i) = Population("pred", 2, 1, a_0P, 1, sigma_gamma, sigma_beta, c_a_pred, c_ss_pred, N0P, p_mut_pred, loci_pred);
        end
    end

    
    
    % Initialize trait values
    for i=1:length(population)
        population(i).trait_values = population(i).update_trait_freq();
        population(i).fitness_values = population(i).update_fitness_values();
    end
    
    % Initialize attack rate for all populations
    prey_trait = [population([population.type]=="prey").individuals.alpha];
    for p = 1:num_populations
        %Prey
        if population(p).type=="prey"
            for i = 1:length(population(p).individuals)
                population(p).individuals(i).a_k = population(p).individuals(i).consumption(resources, num_habitats);
            end
        elseif population(p).type=="pred"
            for i = 1:length(population(p).individuals)
                population(p).individuals(i).g = g;
                population(p).individuals(i).a_k = population(p).individuals(i).consumption(prey_trait, num_habitats);
            end
            population(p).attack_rate = population(p).update_attack_rate(); % Can be redundant and optimized above, to be tested
        end
    end
    
    % Initial variables, for script reproducibility and control
    timestamp = datestr(datetime('now'), 'yymmddHHMMSS');   % Timestamp, used for file naming
    timestamp_day = str2num(extractAfter(timestamp,4));     % Extract day, hour, minute and second as seed
    rng(timestamp_day);                                     % Seed generated with timestamp: DDHHMMSS
   
    % File creation
    if num_populations==2 % When there is only prey
        outfile = strcat('Results/pred_prey_sigmaalpha_',num2str(sigma_alpha), '_pmutprey_',num2str(p_mut_prey), ...
        '_sigmagamma_',num2str(sigma_gamma),'_attack_',num2str(a_0P), ...
        '_g_',num2str(g),'_pmutpred_',num2str(p_mut_pred),'_caprey_',num2str(c_a_prey),'_capred_',num2str(c_a_pred),'_morphsinit_',num2str(morphs),'_sexual_',num2str(is_sexual),'_',timestamp,'.csv');
    else
        outfile = strcat('Results/prey_sigmaalpha_',num2str(sigma_alpha), '_pmutprey_',num2str(p_mut_prey),'_caprey_',num2str(c_a_prey),'_morphsinit_',num2str(morphs),'_sexual_',num2str(is_sexual),'_', timestamp,'.csv');
    end
    outfile_traits = fopen(outfile, 'w');
    

    %% Simulation --------------------------------------------------------
    end_of_simulation = false; % Ends if there's an extinction
    % Iterate for each timestep
    for t = 1:t_end
        
        for i=1:num_populations % For all populations (Prey or pred)
            % Get trait values
            trait_keys = population(i).trait_values.keys();
            trait_val = population(i).trait_values.values();
            fitness_val = population(i).fitness_values.values();
            if length(trait_keys)==0
                end_of_simulation=true;
            end
            % Print trait values and fitness to files
            if mod(t,2)==0
                for j=1:length(trait_keys)
                    if mod(i,2)==1 % Prey file
                        %fprintf('%d\t%.3f\t%d\tprey\n', t, trait_keys{j}, trait_val{j});
                        fprintf(outfile_traits,'%d\t%.3f\t%d\t%d\tprey\n', t, trait_keys{j}, trait_val{j}, fitness_val{j});
                    else
                        %fprintf('%d\t%.3f\t%d\tpred\n', t, trait_keys{j}, trait_val{j});
                        fprintf(outfile_traits,'%d\t%.3f\t%d\t%d\tpred\n', t, trait_keys{j}, trait_val{j}, fitness_val{j});
                    end
                end
            end
        end
       
        if end_of_simulation % If extinction, ends the simulation
            break
        end

        %%% Prey consumption
        prey_pop = [population([population.type]=="prey").individuals];
        cellarr = {prey_pop.a_k};
        a_k_matrix = cat(3,cellarr{:});
        % Calculate equilibrium resource abundance
        for i=1:num_habitats % Get for all current populations
            for j=1:num_resources
                % Calc Rk_eq for the kth resource
                resources(i,j).Rk_eq = resources(i,j).eq_abundance(F, a_k_matrix(i,j,:));
                resource_abundance(t,i,j) = resources(i,j).Rk_eq;
            end
        end

        %%% Predator consumption
        pred_pop = [population([population.type]=="pred").individuals];

        %%% Execute for every population
        for p = 1:num_populations
          % Generate temporal variables for new population
          if population(p).type=="pred"
              new_population = Predator.empty;
          elseif population(p).type=="prey"
              new_population = Prey.empty;
          end

          tmp_popsize = containers.Map('KeyType','double','ValueType','double');
          pop_trait = [population(p).individuals.alpha];

          % Count individuals by trait value
          for i = 1:length(pop_trait)
              if ismember(pop_trait(i), cell2mat(keys(tmp_popsize)))
                tmp_popsize(pop_trait(i)) = tmp_popsize(pop_trait(i)) + 1;
              else
                tmp_popsize(pop_trait(i)) = 1;
              end
          end
          
          %%% Reproduction
          for i = 1:length(population(p).individuals)
              %%% Choose mate, reproduce (or clone) and mutate
              if population(p).type=="prey"
                next_gen = reproduce(i,population(p), F, 0); % NOTE: Prey reproduces asexually, change for sexual reproduction scenario
              else
                next_gen = reproduce(i,population(p), F, is_sexual);
              end

              %%% Offspring
              for j=1:length(next_gen) % For every offspring, disperse (Unused) and evaluate fitness
                  offspring = next_gen(j);
                  if isa(offspring,'Prey') % Prey offspring, calculate new attack rate
                      offspring.a_k = offspring.consumption(resources, num_habitats);
                  end
                  % Disperse (Unused
                  if num_habitats>1 && p_dispersal >= rand()
                      % Switch to a neighbouringh hab
                      offspring.habitat = offspring.disperse(num_habitats);
                      % Recalculate consumption for prey in new habitat
                      if isa(offspring,'Prey')
                        offspring.a_k = offspring.consumption(resources, num_habitats);
                      end
                  end
                  currhab = offspring.habitat;

                  % Calculate fitness
                  pred_attack = zeros(1,i);
                  % Calculate predator attack on individual
                  if width(pred_pop)>0 && isa(offspring,'Prey')
                    pred_attack = [population([population.type]=="pred").attack_rate];
                    pred_attack = pred_attack(:,i);
                  end
                  if num_habitats==1        %% Sympatric model
                      if isa(offspring,'Prey') % Prey offspring
                        offspring.fitness = offspring.calc_fitness_alpha(resources(currhab,:), pred_attack, bmax);
                      else %Predator offspring
                        offspring.a_k = offspring.consumption([prey_pop.alpha], num_habitats);
                        offspring.fitness = offspring.calc_fitness_alpha(bmax);
                      end
                  elseif num_resources==1   %% Allopatric 1-resource model (Not implemented completely)
                    offspring.fitness = offspring.calc_fitness_beta(resources(currhab,:), population.get_popsize(currhab), F);
                  else                      %% Model A+B (Sympatric+allopatric) (Not implemented completely)
                    offspring.fitness = offspring.calc_fitness_full(resources); 
                  end
                  % Survival
                  survival = offspring.fitness/F;
                  if survival>=rand()
                    new_population(length(new_population)+1) = offspring;
                  end
              end
          end
          trait_frequency{t,p} = [cell2mat(keys(tmp_popsize)); cell2mat(values(tmp_popsize))];
          population_size(t,p) = length(new_population);
          population(p).individuals = new_population;
        end
        
        %Update population rait values
        for p=1:num_populations
            population(p).trait_values = population(p).update_trait_freq();
            population(p).fitness_values = population(p).update_fitness_values();
            if population(p).type=="pred"
                prey_trait = [population([population.type]=="prey").individuals.alpha];
                for i = 1:length(population(p).individuals)
                    population(p).individuals(i).a_k = population(p).individuals(i).consumption(prey_trait, num_habitats);
                end
                population(p).attack_rate = population(p).update_attack_rate(); % Can be optimized above, to be tested
            end
        end
        tmp = 1;
        for p=1:num_resources
            for q=1:num_habitats
                tmp = tmp+1;
            end
        end
        for p=1:num_populations
           trait_values = trait_frequency{t,p};
        end
        trait_values = trait_frequency{t,:};
    end
    
    fclose(outfile_traits); % Close file 
end
