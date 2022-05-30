%------------------------------------------------------------
% Julio Ayala
% ju7141ay-s@student.lu.se
% October 2021
% Description: Base class that defines individuals to be part of a
% population. Derives "Predator" and "Prey" objects.
% Usage:
% Create "Individual" objects with
%     Individual(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta, c_a, c_ss, loci, p_mut)
% Where:
%     alpha = alpha trait value
%     beta  = beta trait value
%     a_0 = base attack rate
%     habitat = Number of habitat
%     sigma_alpha = alpha niche width 
%     sigma_beta  = beta niche width 
%     c_a = Alpha trait choosiness parameter
%     c_ss = Sexual trait choosiness parameter
%     loci = Structure with alpha, beta, display and preference genes number of loci
%     p_mut = Mutation rate
%------------------------------------------------------------

classdef Individual < matlab.mixin.Copyable
    properties
        % Traits
        alpha % trait alpha (resource)
        beta % trait beta (habitat)
        dis % Display trait
        pref % Preference trait
        a_0 % base attack rate
        habitat % Curent habitat 
        sigma_alpha % resource niche width
        sigma_beta % habitat niche width
        % Choosiness
        c_a % Choosiness for alpha/gamma value 
        c_ss % Choosines of sexual selection
        % Calculated values
        a_k % attack rate for the kth resource
        fitness % fitness value of individual
        % Genotype
        alpha_gene % Genotype for alpha trait
        beta_gene % Genotype for beta trait
        dis_gene % Genotype for display trait
        pref_gene % Genotype for preference trait
        p_mut % P of mutation
    end
    methods
        %% Constructor
        function obj = Individual(alpha, beta, a_0, habitat, sigma_alpha, sigma_beta, c_a, c_ss, loci, p_mut)
            
            if nargin==0
                obj.alpha = 0;
                obj.beta = 0;
                obj.dis = 0;
                obj.pref = 0;
                obj.a_0 = 0;
                obj.habitat = 0;
                obj.sigma_alpha = 0;
                obj.sigma_beta = 0;
                obj.c_a = 0;
                obj.c_ss = 0;
                obj.p_mut = 0;
            else
                obj.alpha = alpha;
                obj.beta = beta;
                obj.dis = 0;
                obj.pref = 0;
                obj.a_0 = a_0;
                obj.habitat = habitat;
                obj.sigma_alpha = sigma_alpha;
                obj.sigma_beta = sigma_beta;
                obj.c_a = c_a;
                obj.c_ss = c_ss;
                obj.p_mut = p_mut;
            end
            %%% Genotype
            % Ecological traits
            obj.alpha_gene = Genetics('Diallelic', loci.alpha, alpha, 0, 4, p_mut);
            obj.beta_gene = Genetics('Diallelic', loci.beta, beta, 0, 4, p_mut);
            % Display and preference traits
            obj.dis_gene = Genetics('Diallelic', loci.dis, obj.dis, -2, 2, p_mut);
            obj.pref_gene = Genetics('Diallelic', loci.pref, obj.pref, -2, 2, p_mut);
        end
        %% attack rate function
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
        %% Fitness function (alpha/sympatric)
        function fi_alpha = calc_fitness_alpha(obj, resources)
            % Model A
            fi_alpha = sum(obj.a_k.*[resources.Rk_eq]);
        end
        %% Fitness function (beta/allopatric)
        function fi_beta = calc_fitness_beta(obj, resources, nh, F)
            % Model B
            fi_beta = exp(-((obj.beta - obj.habitat)^2)/(2*obj.sigma_beta^2)) * obj.a_0/(1+ ((nh*obj.a_0)/(resources.K.*F)));
        end
        %% Fitness function (Full model)
        function fi = calc_fitness_full(obj, resources)
            rk_eq = zeros(size(resources));
            for i=1:height(resources)
                rk_eq(i,:) = [resources(i,:).Rk_eq];
            end
            %% To fix: Fitness should be based on current habitat
            fi = exp(-((obj.beta - obj.habitat)^2)/(2*obj.sigma_beta^2)) * sum(sum(obj.a_k.*rk_eq));
        end

        %% Dispersal function
        function habitat = disperse(obj, max_hab)
            if obj.habitat==1
                habitat = 2;
            elseif rand()>=0.5 & obj.habitat~=max_hab
                habitat = obj.habitat + 1;
            else
                habitat = obj.habitat - 1;
            end
        end
        
    end
end
