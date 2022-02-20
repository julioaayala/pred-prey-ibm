function main_prey(varargin)
    %Parameters
    t_end = 4000;
    F = 2; % Fecundity rate
    K = 300; % Carrying capacity
    p_dispersal = 0; % Probability of dispersal
    p_mut_prey = 1e-4; % Probability of mutation of prey
    p_mut_pred = 0;
    delta_mut = 0.2; % Mutation delta
    N0N = 100;
    N0P = 10;
    bmax = 1;
    sigma_alpha = 0.5; % Resource niche width
    sigma_gamma = 0.5; % Predator trait niche width
    a_0N = 2; % Base attack rate
    a_0P = 0.008; % Base attack rate
    g = 0.5;
    morphs = 1;

    % Genetics
    loci = 8; % To be implemented

    if ~isempty(varargin)
        sigma_alpha = varargin{1};
        p_mut_prey = varargin{2};
    end
    sigma_beta = 0.5; % Habitat niche width
    num_resources = 3;
    num_habitats = 1;
    num_populations = 1; 
    

    % Initialize resources and populations
    % Resources
    resource_abundance = zeros(t_end,num_habitats, num_resources);
    resources = Resource.empty;
    % Same trait on diff habitats (i)
    for i=1:num_habitats
        for j=1:num_resources
            resources(i,j) = Resource(j,i,K);
        end
    end
    % Individuals
    trait_frequency = {};   
    population_size = zeros(t_end,num_populations);
    fitness = zeros(t_end,num_populations);
    population = Population.empty(num_populations,0);
    % Same resource trait (1) on different habitats with hab trait (i)

    for i=1:num_populations
        if mod(i,2)==1 % Prey
            if morphs>1 % Initialize prey populations with more than 1 trait
                population(i) = Population("prey", 1:morphs, 1, a_0N, 1, sigma_alpha, sigma_beta, N0N, p_mut_prey);
            else
                population(i) = Population("prey", 2, 1, a_0N, 1, sigma_alpha, sigma_beta, N0N, p_mut_prey);
            end
        else
            population(i) = Population("pred", 2, 1, a_0P, 1, sigma_gamma, sigma_beta, N0P, p_mut_pred);
        end
    end
    
    
    % Initialize trait values
    for i=1:length(population)
        population(i).trait_values = population(i).update_trait_freq();
    end
    
    %% Initialize attack rate for the population
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

    
    % Initialize plots
%     figure();
%     legendStr = " Population " + string(1:num_populations);
%     for i=1:num_habitats
%         tmp = "Hab " + i + "; Res " + string(1:num_resources);
%         legendStr = union(legendStr, tmp);
%     end
%     %subplot(2,1,1);
%     colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#80B3FF"];
%     for i=1:(num_populations+(num_resources*num_habitats))
%         h(i) = animatedline("Color", colors(i), "Linewidth",1.5);
%     end
%     legend(legendStr);
%     xlabel("t");
%     ylabel("Abundance");
%     hold on
%     subplot(2,2,3);
%     for i=1:num_populations
%         h2(i) = animatedline("Color", colors(i), 'Marker','*', 'LineStyle','none');
%     end
%     xlabel("t");
%     ylabel("trait value");
% 
%     hold on
%     subplot(2,2,4);
%     h3 = animatedline("Color", colors(5), 'Marker','diamond', 'LineStyle','none');
%     xlabel("Trait");
%     ylabel("Abundance");
%     hold off;

    % Initialize output file
    outfile = strcat('Results/prey_sigmaalpha_',num2str(sigma_alpha), '_', datestr(datetime('now'), 'yymmddHHMMSS'),'.csv');
    outfile_traits = fopen(outfile, 'w');

    %Iterate for each timestep
    % Run simulation
    end_of_simulation = false;
    for t = 1:t_end
        %% Print trait values to files
        if mod(t,2)==0
            for i=1:num_populations
                %Trait values
                trait_keys = population(i).trait_values.keys();
                trait_val = population(i).trait_values.values();
                if length(trait_keys)==0
                    end_of_simulation=true;
                end
                for j=1:length(trait_keys)
                    if mod(i,2)==1 % Prey file
                        %fprintf('%d\t%.3f\t%d\tprey\n', t, trait_keys{j}, trait_val{j});
                        fprintf(outfile_traits,'%d\t%.3f\t%d\tprey\n', t, trait_keys{j}, trait_val{j});
                    else
                        %fprintf('%d\t%.3f\t%d\tpred\n', t, trait_keys{j}, trait_val{j});
                        fprintf(outfile_traits,'%d\t%.3f\t%d\tpred\n', t, trait_keys{j}, trait_val{j});
                    end
                end
            end
        end

        if end_of_simulation
            break
        end

        %% Prey consumption
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

        %% Predator consumption
        pred_pop = [population([population.type]=="pred").individuals];
%         if width(pred_pop)>0
%             cellarr = {pred_pop.a_k};
%             a_k_matrix = cat(3,cellarr{:});
%             % Calculate equilibrium resource abundance
%             for i=1:num_habitats % Get for all current populations
%                 for j=1:num_populations
%                     if population(j).type=="pred"
%                         % Feed
%                     end
%                 end
%             end
%         end

        %% Execute for every population
        for p = 1:num_populations
          % Generate temporal variables for new pop and size
          if population(p).type=="pred"
              new_population = Predator.empty;
          elseif population(p).type=="prey"
              new_population = Prey.empty;
          end
          tmp_popsize = containers.Map('KeyType','double','ValueType','double');
          pop_trait = [population(p).individuals.alpha];
          % Count individuals by habitat
          for i = 1:length(pop_trait)
              if ismember(pop_trait(i), cell2mat(keys(tmp_popsize)))
                tmp_popsize(pop_trait(i)) = tmp_popsize(pop_trait(i)) + 1;
              else
                tmp_popsize(pop_trait(i)) = 1;
              end
          end

          for i = 1:length(population(p).individuals)
              %% Reproduce
              for j=1:F % Create for every offspring
                  offspring = copy(population(p).individuals(i));
                  %% Mutate
                  if population(p).p_mut >= rand()
                    if rand()>0 %%Mutate only resource trait
                        offspring.alpha = offspring.mutate_alpha(delta_mut);
                    else
                        offspring.beta = offspring.mutate_beta(delta_mut);
                    end
                    %% TODO: Fix this
                    if isa(offspring,'Prey') %Prey offspring
                        offspring.a_k = offspring.consumption(resources, num_habitats);
                    end
                  end               
                  %% Disperse
                  if p_dispersal >= rand() && num_habitats>1
                      % Switch to a neighbouringh hab
                      offspring.habitat = offspring.disperse(num_habitats);
                      % Recalculate consumption for prey, since pred is
                      % always calculated
                      if isa(offspring,'Prey') %Prey offspring
                        offspring.a_k = offspring.consumption(resources, num_habitats);
                      end
                  end
                  currhab = offspring.habitat;
                  %% Fitness
                  pred_attack = zeros(i);

                  if width(pred_pop)>0 && isa(offspring,'Prey')
                    pred_attack = [population([population.type]=="pred").attack_rate];
                    pred_attack = pred_attack(:,i);
                  end
                  if num_habitats==1 % Model A
                      if isa(offspring,'Prey') %Prey offspring
                        offspring.fitness = offspring.calc_fitness_alpha(resources(currhab,:), pred_attack, bmax);
                      else %Predator offspring
                        offspring.a_k = offspring.consumption([prey_pop.alpha], num_habitats);
                        offspring.fitness = offspring.calc_fitness_alpha(bmax);
                      end
                  elseif num_resources==1 % Model B
                    offspring.fitness = offspring.calc_fitness_beta(resources(currhab,:), population.get_popsize(currhab), F);
                  else % Model A+B
                    offspring.fitness = offspring.calc_fitness_full(resources); 
                  end
                  survival = offspring.fitness/F;
                  %% Survive
                  if survival>=rand()
                    new_population(length(new_population)+1) = offspring;
                  end
              end
          end
          trait_frequency{t,p} = [cell2mat(keys(tmp_popsize)); cell2mat(values(tmp_popsize))];
          population_size(t,p) = length(new_population);
          population(p).individuals = new_population;
        end
        
        %Update figures and trait values
        for p=1:num_populations
            population(p).trait_values = population(p).update_trait_freq();
            if population(p).type=="pred"
                %if isempty(population(p).individuals)
                %    disp("RIP");
                %end
                prey_trait = [population([population.type]=="prey").individuals.alpha];
                for i = 1:length(population(p).individuals)
                    population(p).individuals(i).a_k = population(p).individuals(i).consumption(prey_trait, num_habitats);
                end
                population(p).attack_rate = population(p).update_attack_rate(); % Can be redundant and optimized above, to be tested
            end
%             addpoints(h(p),t,population_size(t,p));
        end
        % Add to figure p+num_populations
        tmp = 1;
        for p=1:num_resources
            for q=1:num_habitats
%                 addpoints(h(tmp+num_populations),t,resource_abundance(t,q,p));
                tmp = tmp+1;
            end
        end
        for p=1:num_populations
           trait_values = trait_frequency{t,p};
%            addpoints(h2(p),repelem(t,length(trait_values)),trait_values);
        end

        trait_values = trait_frequency{t,:};
        %clearpoints(h3);
        %addpoints(h3,trait_values(1,:),trait_values(2,:));
        %addpoints(h3,trait_values(1,trait_values(2,:)>50),trait_values(2,trait_values(2,:)>50));
        %drawnow;
    end
    %fclose(preyfile_traits);
    fclose(outfile_traits);
    %exit();
end
