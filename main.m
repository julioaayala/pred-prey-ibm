function main(varargin)
    %Parameters
    t_end = 1000;
    F = 2; % Fecundity rate
    K = 300; % Carrying capacity
    P_dispersal = 0.2; % Probability of dispersal
    P_mutate = 0; % Probability of dispersal
    delta_mut = 0.4; % Mutation delta
    a_0 = 2; % Base attack rate
    sigma_alpha = 0.7; % Resource niche width
    sigma_beta = 0.7; % Habitat niche width
    num_resources = 2;
    num_habitats = 1;
    num_populations = 2;

    % Initialize resources and populations
    % Resources
    resource_abundance = zeros(t_end,num_habitats, num_resources);
    resources = Resource.empty;
    % Same trait on diff habitats (i)
    for i=1:num_habitats
        for j=1:num_resources
            % habitat = i; trait = j; K
            if j==2
                resources(i,j) = Resource(j,i,K);
            
            else
                resources(i,j) = Resource(j,i,K);
            end
        end
    end
    % Individuals
    trait_frequency = {};   
    population_size = zeros(t_end,num_populations);
    fitness = zeros(t_end,num_populations);
    population = Population.empty(num_populations,0);
    % Same resource trait (1) on different habitats with hab trait (i)
%     for p=1:num_populations
%         % alpha = p, beta = i a_0 = a_0, habitat = i, 
%         % sigma_alpha = sigma_alpha
%         % Initialize ech population with 2 individuals
%         if mod(p,2)==0
%             population(p) = Population("pred", 1, 1, a_0, 1, sigma_alpha, sigma_beta, 2);
%         else
%             population(p) = Population("prey", 1, 1, a_0, 1, sigma_alpha, sigma_beta, 4);
%         end
%     end

    population(1) = Population("prey", 1, 1, a_0, 1, sigma_alpha, sigma_beta, 2);
    population(2) = Population("pred", 1, 1, a_0/5, 1, sigma_alpha, sigma_beta, 2);
    population(3) = Population("prey", 2, 1, a_0, 1, sigma_alpha, sigma_beta, 2);
    
    %% Initialize attack rate for the population
    for p = 1:num_populations
        %Prey
        if population(p).type=="prey"
            for i = 1:length(population(p).individuals)
                population(p).individuals(i).a_k = population(p).individuals(i).consumption(resources, num_habitats);
            end
        elseif population(p).type=="pred"
            prey_abundances = population([population.type]=="prey").get_abundances("alpha");
            for i = 1:length(population(p).individuals)
                population(p).individuals(i).a_k = population(p).individuals(i).consumption(prey_abundances);
            end
        end
    end

    
    % Initialize plots
    figure();
    legendStr = " Population " + string(1:num_populations);
    for i=1:num_habitats
        tmp = "Hab " + i + "; Res " + string(1:num_resources);
        legendStr = union(legendStr, tmp);
    end
    subplot(2,1,1);
    colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#80B3FF"];
    for i=1:(num_populations+(num_resources*num_habitats))
        h(i) = animatedline("Color", colors(i));
    end
    legend(legendStr);
    xlabel("t");
    ylabel("abundance");
    hold on
    subplot(2,2,3);
    for i=1:num_populations
        h2(i) = animatedline("Color", colors(i), 'Marker','*', 'LineStyle','none');
    end
    xlabel("t");
    ylabel("trait value");

    hold on
    subplot(2,2,4);
    h3 = animatedline("Color", colors(5), 'Marker','diamond', 'LineStyle','none');
    xlabel("Trait");
    ylabel("Abundance");
    hold off;

    %Iterate for each timestep
    % Run simulation
    for t = 1:t_end
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
        if width(pred_pop)>0
            cellarr = {pred_pop.a_k};
            a_k_matrix = cat(3,cellarr{:});
            % Calculate equilibrium resource abundance
            for i=1:num_habitats % Get for all current populations
                for j=1:num_populations
                    if population(j).type=="pred"
                        % Feed
                    end
                end
            end
        end

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
                  if P_mutate >= rand()
                    if rand()>0.5
                        offspring.alpha = offspring.mutate_alpha(delta_mut);
                    else
                        offspring.beta = offspring.mutate_beta(delta_mut);
                    end
                    %% TODO: Fix this
                    offspring.a_k = offspring.consumption(resources, num_habitats);
                  end               
                  %% Disperse
                  if P_dispersal >= rand() & num_habitats>1
                      % Switch to a neighbouringh hab
                      offspring.habitat = offspring.disperse(num_habitats);
                      % Recalculate consumption
                      offspring.a_k = offspring.consumption(resources, num_habitats);
                  end
                  currhab = offspring.habitat;
                  %% Fitness
                  prey_abundances = population([population.type]=="prey").get_abundances("alpha");
                  pred_abundances = 0;
                  pred_attack = 0;

                  if width(pred_pop)>0
                    new_prey_abundances = population([population.type]=="pred").get_eq_abundance_prey(prey_abundances, F);
                    pred_abundances = population([population.type]=="pred").get_abundances("alpha");
                    pred_attack = population([population.type]=="pred").sum_ak();
                  end
                  if num_habitats==1 % Model A
                      if isa(offspring,'Prey') %Prey offspring
                        offspring.fitness = offspring.calc_fitness_alpha(resources(currhab,:),prey_abundances,pred_abundances,pred_attack);
                      else %Predator offspring
                          disp(new_prey_abundances);
                          offspring.fitness = offspring.calc_fitness_alpha(new_prey_abundances);
                      end
                  elseif num_resources==1 % Model B
                    offspring.fitness = offspring.calc_fitness_beta(resources(currhab,:), population.get_popsize(currhab), F);
                  else % Model A+B
                    offspring.fitness = offspring.calc_fitness_full(resources); 
                  end
                  survival = offspring.fitness/F;
                  %% Survive
                  if survival>=rand()
                      if isa(offspring,'Prey')
                        % If prey is not eaten (i.e. ratio between new
                        % abundance and previous abundance)
                        if width(pred_pop)>0
                            p_noteaten = new_prey_abundances(offspring.habitat, offspring.alpha)/prey_abundances(offspring.habitat, offspring.alpha);
                        else
                            p_noteaten=1;
                        end
                        if p_noteaten>rand()
                            new_population(length(new_population)+1) = offspring;
                        end
                      else
                        new_population(length(new_population)+1) = offspring;
                      end
                  end
              end
          end
          trait_frequency{t,p} = [cell2mat(keys(tmp_popsize)); cell2mat(values(tmp_popsize))];
          population_size(t,p) = length(new_population);
          population(p).individuals = new_population;
          if size(new_population)==0
             disp("RIP");
          end
        end
        
        %Update figures
        for p=1:num_populations
           addpoints(h(p),t,population_size(t,p));
        end
        % Add to figure p+num_populations
        tmp = 1;
        for p=1:num_resources
            for q=1:num_habitats
                addpoints(h(tmp+num_populations),t,resource_abundance(t,q,p));
                tmp = tmp+1;
            end
        end
        for p=1:num_populations
           trait_values = trait_frequency{t,p};
           addpoints(h2(p),repelem(t,length(trait_values)),trait_values);
        end

        trait_values = trait_frequency{t,:};
        clearpoints(h3);
        %addpoints(h3,trait_values(1,:),trait_values(2,:));
        addpoints(h3,trait_values(1,trait_values(2,:)>50),trait_values(2,trait_values(2,:)>50));
        drawnow;
    end


end