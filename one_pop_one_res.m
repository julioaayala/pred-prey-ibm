function one_pop_one_res(varargin)
    
    %Parameters
    t_end = 200;
    F = 2; % Fecundity rate
    K = 300; % Carrying capacity
    P_d = 0; % Probability of dispersal
    a_0 = 2; % Base attack rate
    sigma_alpha = 1; % Resource niche width
    sigma_beta = 1; % Habitat niche width
    num_resources = 1;
    num_habitats = 1;
    num_populations = 1;

    % Initialize resources and populations
    % Resources
    resource_abundance = zeros(t_end,num_habitats, num_resources);
    resources = Resource.empty;
    % Same trait on diff habitats (i)
    for i=1:num_habitats
        for j=1:num_resources
            % habitat = i; trait = j; K
            resources(i,j) = Resource(j,i,K);
        end
    end
    % Individuals
    population_size = zeros(t_end,num_populations);
    population = Individual.empty(num_populations,0);
    % Same resource trait (1) on different habitats with hab trait (i)
    for i=1:num_populations
        % alpha = 1, beta = i a_0 = a_0, habitat = i, 
        % sigma_alpha = sigma_alpha
        population(i) = Individual(i, 1, a_0, 1, sigma_alpha, sigma_beta);
    end
    % Initialize attack rate for the population
    for i = 1: length(population)
        population(i).a_k = population(i).consumption(resources, num_habitats);
    end
    %Iterate for each timestep
    
    % Initialize plot
    figure;
    legendStr = union("Population " + string(1:num_populations), ...
                      "Resource " + string(1:num_resources));
    
    colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"];
    for i=1:(num_populations+num_resources)
        h(i) = animatedline("Color", colors(i));
    end
    legend(legendStr);
    xlabel("t");
    ylabel("abundance");

    % Run simulation
    for t = 1:t_end
        % Get attack rate from population in matrix form
        cellarr = {population.a_k};
        
        a_k_matrix = cat(3,cellarr{:});

        % Calculate equilibrium resource abundance
        
        for i=1:max([population.habitat]) % Get for all current populations
            for j=1:num_resources
                % Calc Rk_eq for the kth resource
                resources(i,j).Rk_eq = resources(i,j).eq_abundance(F, a_k_matrix(i,j,:));
                resource_abundance(t,i,j) = resources(i,j).Rk_eq;
            end
        end
        
        % Generate temporal variables for new pop and size
        new_population = Individual.empty;
        tmp_popsize = zeros(1,num_populations);
        pop_trait = [population.alpha];
        % Count individuals by habitat
        for i = 1:length(pop_trait)
            tmp_popsize(pop_trait(i)) = tmp_popsize(pop_trait(i)) + 1;
        end
        for i = 1:length(population)
            % Reproduction
            for j=1:F % Create for every offspring
                offspring = copy(population(i));
                % Dispersal
                if P_d >= rand()
                    % Switch to a neighbouringh hab
                    offspring.habitat = offspring.disperse(num_habitats);
                    % Recalculate consumption
                    offspring.a_k = offspring.consumption(resources, num_habitats);
                end
                currhab = offspring.habitat;
                offspring.fitness = offspring.calc_fitness_alpha(resources(currhab,:));
                survival = offspring.fitness/F;
                if survival>=rand()
                    new_population(length(new_population)+1) = offspring;
                end
            end
        end
        population_size(t,:) = tmp_popsize;
        population = new_population;
        
        %Update figures
        for p=1:num_populations
            addpoints(h(p),t,population_size(t,p));
        end
        % Add to figure p+num_populations
        for p=1:num_resources
            addpoints(h(p+num_populations),t,resource_abundance(t,1,p));
        end

        drawnow;
    end


end