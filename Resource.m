%------------------------------------------------------------
% Julio Ayala
% ju7141ay-s@student.lu.se
% October 2021
% Description: Resource class for ecological interactions
% Usage:
% Create resource objects from a script with
%     Resource(trait, habitat, K)
% Where:
%     trait = Ecological trait value
%     habitat = Habitat trait value
%     K = System size
%------------------------------------------------------------

classdef Resource
    properties
        trait       % k Trait [1..3]
        habitat     % Habitat number #
        K           % Carrying capacity
        Rk_eq       % Equilibrium resouce abundance 
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
        % Function to calculate equilibrium resource abundance from
        % consumer attack
        function Rk_eq = eq_abundance(obj, F, a_ik) 
            Rk_eq = obj.K / (1+sum(a_ik./F));
        end
    end
end


