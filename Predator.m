classdef Predator < Individual
    properties
        g; % Conversion factor
    end
    methods
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
        function fi_alpha = calc_fitness_alpha(obj, bmax)
            % Model A
            fi_alpha = obj.g * bmax * sum(obj.a_k);
        end
    end
end