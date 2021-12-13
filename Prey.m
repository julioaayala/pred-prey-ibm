classdef Prey < Individual
    methods
        function fi_alpha = calc_fitness_alpha(obj, resources, pred_attack, bmax)
            % Model A
            a_k = obj.a_k(obj.habitat,:);
            fi_alpha = sum(a_k.*[resources.Rk_eq]) - bmax * sum(pred_attack);
        end
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
    end
end