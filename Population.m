classdef Population
    properties
        individuals
        type
        aik
    end
    methods
        % Constructor
        function obj = Population(type, alpha, beta, a_0, habitat, sigma_alpha, sigma_beta, n)
            if nargin==0
                obj.individuals = Individual.empty(1,0);
            else
                obj.type = type;
                if type=="pred"
                    obj.individuals = Predator.empty(n,0);
                elseif type=="prey"
                    obj.individuals = Prey.empty(n,0);
                end
                for i=1:n
                    if type=="pred"
                        obj.individuals(i) = Predator(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta);
                    elseif type=="prey"
                        obj.individuals(i) = Prey(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta);
                    end
                end
            end
        end
        % Get the frequency of a habitat in populations.
        function pop_size_habitat = get_popsize(obj, habitat)
            habs = [obj.individuals.habitat];
            pop_size_habitat = sum(habs(:)==habitat);
        end
        function abundances = get_abundances(obj, trait)
            pop = [obj.individuals];
            disp(pop);
            hab_list = [pop.habitat];
            if trait=="alpha"
                trait_list = [pop.alpha];
            elseif trait=="beta"
                trait_list = [pop.beta];
            end 
            abundances = zeros(max(hab_list),max(trait_list));
            disp(trait_list);
            for i=1:max(hab_list)
                for j=1:max(trait_list)
                    abundances(i,j) = nnz(hab_list==i & trait_list==j);
                end
            end
        end
        function eq_abundances = get_eq_abundance_prey(obj, abundances, F)
            eq_abundances = zeros(size(abundances));
            cellarr = {obj.individuals.a_k};
            a_k_matrix = cat(3,cellarr{:});
            for i=1:height(abundances)
                for j=1:height(abundances)
                    eq_abundances(i,j) = abundances(i,j)/(1+sum(a_k_matrix(i,j,:)/F));
                end
            end
        end

        function sum_ak = sum_ak(obj)
            prey_pop = [obj.individuals];
            cellarr = {prey_pop.a_k};
            a_k_matrix = cat(3,cellarr{:});
            sum_ak = sum(a_k_matrix,3);
        end
    end
end
