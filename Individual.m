classdef Individual < matlab.mixin.Copyable
    properties
        alpha % trait alpha (resource)
        beta % trait beta (habitat)
        a_0 % base attack rate
        habitat % Curent habitat 
        sigma_alpha % resource niche width
        sigma_beta % habitat niche width
        a_k % attack rate for the kth resource
        fitness % fitness value of individual
        alpha_genotype % Genotype for alpha trait
    end
    methods
        % Constructor
        function obj = Individual(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta)
            if nargin==0
                obj.alpha = 0;
                obj.beta = 0;
                obj.a_0 = 0;
                obj.habitat = 0;
                obj.sigma_alpha = 0;
                obj.sigma_beta = 0;
            else
                obj.alpha = alpha;
                obj.beta = beta;
                obj.a_0 = a_0;
                obj.habitat = habitat;
                obj.sigma_alpha = sigma_alpha;
                obj.sigma_beta = sigma_beta;
                obj.alpha_genotype = Genetics('Diallelic', 8, alpha, 0, 4);
            end
        end
        % attack rate function
        function aik = consumption(obj, resources, habitats)
            % a/K -> max attack rate per resource density unit
            % a consumer i specialized on resource k has a_i = k
            % sigma -> niche width
            aik = zeros(habitats,max([resources.trait]));
            for j=1:height(resources)
                % Only assign consumption to habitat the ind. is in
                for k=1:length(resources(j,:))
                    if resources(j,k).habitat == obj.habitat
                        aik(obj.habitat,resources(j,k).trait) = (obj.a_0/resources(j,k).K)*exp(-((obj.alpha-resources(j,k).trait)^2)/(2*obj.sigma_alpha^2));
                    end
                end
            end
        end
        % Fitness function (alpha)
        function fi_alpha = calc_fitness_alpha(obj, resources)
            % Model A
            fi_alpha = sum(obj.a_k.*[resources.Rk_eq]);
        end
        % Fitness function (beta)
        function fi_beta = calc_fitness_beta(obj, resources, nh, F)
            % Model B
            fi_beta = exp(-((obj.beta - obj.habitat)^2)/(2*obj.sigma_beta^2)) * obj.a_0/(1+ ((nh*obj.a_0)/(resources.K.*F)));
        end
        % Fitness function (Full model)
        function fi = calc_fitness_full(obj, resources)
            rk_eq = zeros(size(resources));
            for i=1:height(resources)
                rk_eq(i,:) = [resources(i,:).Rk_eq];
            end
            %% TODO Should fitness be based on curr hab only?
            fi = exp(-((obj.beta - obj.habitat)^2)/(2*obj.sigma_beta^2)) * sum(sum(obj.a_k.*rk_eq));
        end
        % Dispersal function
        function habitat = disperse(obj, max_hab)
            if obj.habitat==1
                habitat = 2;
            elseif rand()>=0.5 & obj.habitat~=max_hab
                habitat = obj.habitat + 1;
            else
                habitat = obj.habitat - 1;
            end
        end
        function subpop = getsubpopalpha(obj, trait)
            subpop = Individual.empty;
            for i=1:length(obj)
                if obj(i).alpha == trait
                    subpop(length(subpop)+1) = obj(i);
                end
            end
        end
        function new_trait = mutate_alpha(obj, delta)
            if rand() >= 0.5
                new_trait = obj.alpha+delta; 
            else 
                new_trait = obj.alpha-delta; 
            end
        end
        function new_trait = mutate_beta(obj, delta)
            if rand() >= 0.5
                new_trait = obj.beta+delta; 
            else 
                new_trait = obj.beta-delta; 
            end
        end
    end
end
