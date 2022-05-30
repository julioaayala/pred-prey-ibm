%------------------------------------------------------------
% Julio Ayala
% ju7141ay-s@student.lu.se
% December 2021
% Description: Genetics class that codes for one diallelic continuous trait
% Usage:
% Create gene objects with
%     Genetics(type, loci, value, minval, maxval, p_mut)
% Where:
%     type = "diallelic" (Only option available now)
%     loci = Number of loci
%     value = Phenotype trait value
%     minval, maxval = Range of values for the trait
%     p_mut = Mutation rate
%------------------------------------------------------------

classdef Genetics < matlab.mixin.Copyable
    properties
        type    % Diallelic, for now
        loci    % Number of loci
        genotype % Matrix of binary values
        minval  % range of phenotype values
        maxval
        p_mut
    end

    methods
        %% Constructor
        function obj = Genetics(type, loci, value, minval, maxval, p_mut)
            if nargin==0
                obj.type = "Diallelic";
                obj.loci = 1;
            else
                obj.type = type;
                obj.loci = loci;
                obj.minval = minval;
                obj.maxval = maxval;
                obj.p_mut = p_mut;
            end
            obj.genotype = phen_to_gen(obj, value);
        end
        
        %% Phenotype to genotype function, to be used only for initialization
        function genotype = phen_to_gen(obj, value)
            if obj.type == "Diallelic"
                genotype = zeros([obj.loci,2]);
                % Rescale phen value
                rescaled = obj.loci*2*(value - obj.minval)/(obj.maxval - obj.minval) + obj.loci*-1;
                temp_phen = obj.loci*-1;
                i = 1;
                j = 0;
                % Modify genotype until it reaches the rescaled value
                while temp_phen<rescaled
                    genotype(i,mod(j,2)+1)=1;
                    if mod(j,2)==1
                        i = i+1;
                    end
                    j = j+1;
                    temp_phen = temp_phen+1;
                end
            end
        end
        
        %% Genotype to phenotype function
        function phenotype = gen_to_phen(obj)
            % Calculate phenotype based on 
            rescaled = sum(sum(obj.genotype==1))/2 - sum(sum(obj.genotype==0))/2;
            % Normalizing 
            phenotype = obj.minval + (obj.maxval-obj.minval) * (rescaled-(obj.loci*-1))/(2*obj.loci);
        end
        
        %% Mutation function 
        function obj = mutate(obj)
            % Function to mutate on each loci based on p_mut
            rangen = rand(obj.loci,1);
            mutations = find(obj.p_mut>=rangen); % Get indices of mutations
            for i=1:length(mutations)
                j = randsample(2,1); % allele
                obj.genotype(mutations(i),j) = abs(obj.genotype(mutations(i),j)-1); % Mutate locus
            end
        end
        
        %% Recombination function
        function obj = recombinate(obj, obj2)
            recombined = zeros([obj.loci,2]);
            % Choose 1 allele from each parent
            for i=1:obj.loci
                recombined(i,:) = [obj.genotype(i,randsample(2,1)), obj2.genotype(i,randsample(2,1))];
            end
            obj.genotype = recombined;
        end
    end
end
