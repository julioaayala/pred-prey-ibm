%------------------------------------------------------------
% Julio Ayala
% ju7141ay-s@student.lu.se
% December 2021
% Description: Predator class to model prey individuals, inherits from Individual
% Contains functions for resource consumption and fitness
%------------------------------------------------------------

classdef Predator < Individual
    properties
        g; % Conversion factor
    end
    methods
        %% Function to calculate resource consumption (i.e. attack rate)
        function aik = consumption(obj, prey_trait, num_habitats)
            % a/K -> max attack rate per resource density unit
            % a consumer i specialized on resource/prey k has a_i = k
            % sigma -> niche width
            aik = zeros(num_habitats, length(prey_trait));
            for h=1:num_habitats
                if h==obj.habitat
                    for i=1:length(prey_trait)
                        aik(h,i) = obj.a_0*exp(-((obj.alpha-prey_trait(i))^2)/(2*obj.sigma_alpha^2));
                    end
                end
            end
        end
        %% Function to calculate fitness in sympatry
        function fi_alpha = calc_fitness_alpha(obj, bmax)
            % Model A
            fi_alpha = obj.g * bmax * sum(obj.a_k);
        end
    end
end