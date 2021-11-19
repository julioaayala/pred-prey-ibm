classdef Resource
    properties
        trait % k Trait [1..5]
        habitat % Hab #
        K % Carrying capacity
        Rk_eq % Eq. resouce abundance 
    end
    methods
        % Constructor
        function obj = Resource(trait, habitat, K)
            if nargin==0
                obj.trait = 0;
                obj.habitat = 0;
                obj.K = 0;
            else
                obj.trait = trait;
                obj.habitat = habitat;
                obj.K = K;
            end
        end
        % Function for equilibrium resource abundance
        function Rk_eq = eq_abundance(obj, F, a_ik) 
            Rk_eq = obj.K / (1+sum(a_ik./F));
        end
    end
end


