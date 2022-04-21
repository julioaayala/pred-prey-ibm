function next_gen = reproduce(ind_index, pop, F, is_sexual)
    % Function to reproduce, generating F offspring
    ind = pop.individuals(ind_index);
    % Prepare offspring array
    if pop.type=="pred"
        next_gen = Predator.empty;
    elseif pop.type=="prey"
        next_gen = Prey.empty;
    end
    if is_sexual 
        %% Sexual case
        %% Mate choice
        % Choose from up to 100 individuals (Or pop size)
        maxmates = min(1000,length(pop.individuals)-1);
        % Exclude the individual reproducing
        popexcl = pop.individuals(setdiff(1:length(pop.individuals),ind_index));
        % Sample possible mates without replacement
        samples = randsample(popexcl, maxmates, false);
        i = 1;
        found_mate = false;
        % Browse in pool while it doesn't find a mate
        while ~(found_mate || i>maxmates)
            ind2 = samples(i);
            % Probability of acceptance by trait alpha
            P_alpha = p_assortative(ind.c_a, ind.alpha, ind2.alpha);
            % Probability of acceptance by preference/display traits *TBI
            % P_ss = p_assortative(ind.c_ss, ind.pref, ind2.dis);
            P_ss = 1;
            % Joint probability
            P_A = P_alpha * P_ss;
            if P_A>= rand()
                found_mate = true;
            else
                i = i+1;
            end
        end
    
        % Reproduce only if it finds a mate
        if found_mate
            for i=1:F
                offspring = copy(ind);
                offspring.alpha_gene = copy(ind.alpha_gene);
                offspring.dis_gene = copy(ind.dis_gene);
                offspring.pref_gene = copy(ind.pref_gene);
                % Recombine
                offspring.alpha_gene.recombinate(ind2.alpha_gene);
                offspring.dis_gene.recombinate(ind2.dis_gene);
                offspring.pref_gene.recombinate(ind2.pref_gene);
                % Probably mutate
                offspring.alpha_gene.mutate();
                offspring.dis_gene.mutate();
                offspring.pref_gene.mutate();
                % Get phenotype
                offspring.alpha = offspring.alpha_gene.gen_to_phen();
                offspring.dis = offspring.dis_gene.gen_to_phen();
                offspring.pref = offspring.pref_gene.gen_to_phen();
                next_gen(i) = offspring;
            end
        end
    else 
        %% Asexual case
        for i=1:F
            offspring = copy(ind);
            %% Probably mutate
            offspring.alpha_gene = copy(ind.alpha_gene);
            offspring.dis_gene = copy(ind.dis_gene);
            offspring.pref_gene = copy(ind.pref_gene);
            offspring.alpha_gene.mutate();
            offspring.dis_gene.mutate();
            offspring.pref_gene.mutate();
            % Get phenotype
            offspring.alpha = offspring.alpha_gene.gen_to_phen();
            offspring.dis = offspring.dis_gene.gen_to_phen();
            offspring.pref = offspring.pref_gene.gen_to_phen();
            next_gen(i) = offspring;
        end
    end
end

% Function to evaluate compatibility between two individuals 
function p_accept = p_assortative(c_a, trait1, trait2)
    if c_a>=0
        p_accept = exp(-c_a*((trait2-trait1)^2));
    else
        p_accept = min(1,exp(-c_a*(((trait2-trait1)^2)-1)));
    end
end