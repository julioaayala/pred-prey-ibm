classdef Genetics < matlab.mixin.Copyable
    properties
        type    % Diallelic, for now
        loci    % Number of loci
        genotype % Matrix of binary values
        minval  % range of phenotype values
        maxval  
    end
    methods
        % Constructor
        function obj = Genetics(type, loci, value, minval, maxval)
            if nargin==0
                obj.type = "Diallelic";
                obj.loci = 1;
            else
                obj.type = type;
                obj.loci = loci;
                obj.minval = minval;
                obj.maxval = maxval;
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
        %% Mutate
        function obj = mutate_genotype(obj, p_mut)
            if p_mut >= rand()
                i = randsample(obj.loci,1); % locus pos
                j = randsample(2,1); % allele pos
                obj.genotype(i,j) = abs(obj.genotype(i,j)-1); % Mutate locus
            end
        end
        %% Recombination
        function recombined = recombinate(obj, obj2)
            i = randsample(obj.loci,1); % locus split position
            j = randsample(2,1); % allele pos
        end
    end
end
