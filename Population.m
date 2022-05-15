classdef Population
    properties
        individuals
        type
        trait_values
        fitness_values
        attack_rate
    end
    methods
        % Constructor
        function obj = Population(type, alpha, beta, a_0, habitat, sigma_alpha, sigma_beta, c_a, c_ss, n, p_mut, loci)
            if nargin==0
                obj.individuals = Individual.empty(1,0);
            else
                obj.type = type;
                if type=="pred"
                    obj.individuals = Predator.empty(n*length(alpha),0);
                elseif type=="prey"
                    obj.individuals = Prey.empty(n*length(alpha),0);
                end
                %% Polimorphic initial population
                if length(alpha)>1
                    for a=1:length(alpha)
                        for i=1:n
                            if type=="pred"
                                obj.individuals((a-1)*n+i) = Predator(alpha(a), beta, a_0, habitat, sigma_alpha, sigma_beta, c_a, c_ss, loci, p_mut);
                            elseif type=="prey"
                                obj.individuals((a-1)*n+i) = Prey(alpha(a), beta, a_0, habitat, sigma_alpha, sigma_beta, c_a, c_ss, loci, p_mut);
                            end
                        end
                    end
                %% Monomorphic initical population
                else
                    for i=1:n
                        if type=="pred"
                            obj.individuals(i) = Predator(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta, c_a, c_ss, loci, p_mut);
                        elseif type=="prey"
                            obj.individuals(i) = Prey(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta, c_a, c_ss, loci, p_mut);
                        end
                    end
                end
            end
        end
        % Get the frequency of a habitat in populations.
        function pop_size_habitat = get_popsize(obj, habitat)
            habs = [obj.individuals.habitat];
            pop_size_habitat = sum(habs(:)==habitat);
        end
        % Get the frequency of a trait in populations.
        function abundances = get_abundances(obj, trait)
            pop = [obj.individuals];
            hab_list = [pop.habitat];
            if trait=="alpha"
                trait_list = [pop.alpha];
            elseif trait=="beta"
                trait_list = [pop.beta];
            end 
            abundances = zeros(max(hab_list),max(trait_list));
            for i=1:max(hab_list)
                for j=1:max(trait_list)
                    abundances(i,j) = nnz(hab_list==i & trait_list==j);
                end
            end
        end
        % Get the total attack rate of the population
        function sum_ak = sum_ak(obj)
            prey_pop = [obj.individuals];
            cellarr = {prey_pop.a_k};
            a_k_matrix = cat(3,cellarr{:});
            sum_ak = sum(a_k_matrix,3);
        end
        % Function to update the trait frequency attribute of the
        % population
        function trait_freq = update_trait_freq(obj)
            alpha = [obj.individuals.alpha];
            unique_vals = unique(alpha);
            trait_freq = containers.Map('KeyType','double','ValueType','double');
            for i=1:length(unique_vals)
                trait_freq(unique_vals(i)) = length(find(alpha==unique_vals(i)));
            end
        end
        % Function to update the fitness attribute of the
        % population
        function trait_fit = update_fitness_values(obj)
            alpha = [obj.individuals.alpha];
            unique_vals = unique(alpha);
            trait_fit = containers.Map('KeyType','double','ValueType','double');
            for i=1:length(unique_vals)
                subpop = obj.individuals([obj.individuals.alpha]==unique_vals(i));
                trait_fit(unique_vals(i)) = mean([subpop.fitness]);
            end
        end
        % Function to update the attack rate attribute of the
        % population
        function attack_rate = update_attack_rate(obj)
          if ~isempty(obj.individuals)
            attack_rate = zeros(length(obj.individuals), length(obj.individuals(1).a_k));
            for i=1:length(obj.individuals)
              attack_rate(i,:) = obj.individuals(i).a_k;
            end
          else
            attack_rate = zeros(1,1);
          end
        end
    end
end
