classdef Population
    properties
        individuals
        aik
    end
    methods
        % Constructor
        function obj = Population(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta, n)
            if nargin==0
                obj.individuals = Individual.empty(1,0);
            else
                obj.individuals = Individual.empty(n,0);
                for i=1:n
                    obj.individuals(i) = Individual(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta);
                end
            end
        end
        % Get the frequency of a habitat in populations.
        function pop_size_habitat = get_popsize(obj, habitat)
            habs = [obj.individuals.habitat];
            pop_size_habitat = sum(habs(:)==habitat);
        end
    end
end
