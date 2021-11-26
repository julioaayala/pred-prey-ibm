classdef Prey < Individual
    methods
        function fi_alpha = calc_fitness_alpha(obj, resources, prey_pop, pred_pop,pred_attack)
            % Model A
            fi_alpha = sum(obj.a_k.*[resources.Rk_eq]) - sum(pred_pop(obj.habitat,:)./prey_pop(obj.habitat,:));
        end
    end
end