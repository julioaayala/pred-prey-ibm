classdef Predator < Individual
    properties
        g = 0.4; % Conversion factor
    end
    methods
        function aik = consumption(obj, prey_abundance)
            % a/K -> max attack rate per resource density unit
            % a consumer i specialized on resource/prey k has a_i = k
            % sigma -> niche width
            aik = zeros(height(prey_abundance),width(prey_abundance));
            for i=1:height(prey_abundance)
                if i==obj.habitat 
                    % Only assign consumption to habitat the ind. is in
                    for j=1:length(prey_abundance(i,:))
                        aik(i,j) = (obj.a_0/prey_abundance(i,j))*exp(-((obj.alpha-j)^2)/(2*obj.sigma_alpha^2));
                    end
                end
            end
        end
        function fi_alpha = calc_fitness_alpha(obj, prey_abundances)
            % Model A
            fi_alpha = obj.g * sum(obj.a_k.*prey_abundances(obj.habitat,:));
        end
    end
end