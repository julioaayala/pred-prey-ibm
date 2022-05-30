function fitness_landscape_test(varargin)
    alpha = 1:3;
    gamma = 1.5:0.0625:2.50;
    sigma_gamma = 0.3:0.1:0.6;
    resources = [Resource(1,1,400), Resource(2,1,400), Resource(3,1,400)];
    loci_prey.alpha = 16;
    loci_prey.beta = 8;
    loci_prey.dis = 8;
    loci_prey.pref = 8;

    loci_pred.alpha = 32;
    loci_pred.beta = 8;
    loci_pred.dis = 8;
    loci_pred.pref = 8;

    preypop = Population("prey", 2, 1, 2, 1, 0.45, 1, 1, 0, length(alpha)*400, 1e-5, loci_prey);
    predpop = Population("pred", 2, 1, 0.005, 1, 0.2, 1, 1, 1, length(gamma)*length(sigma_gamma), 1e-4, loci_pred);

    for i=1:length(alpha)
      for j=1:100  
        preypop.individuals(((i-1)*100)+j).alpha = alpha(i);
      end
    end

    for i=1:length(gamma)
        for j=1:length(sigma_gamma)
            predpop.individuals(((i-1)*length(sigma_gamma))+j).g = 0.65;
            predpop.individuals(((i-1)*length(sigma_gamma))+j).alpha = gamma(i);
            predpop.individuals(((i-1)*length(sigma_gamma))+j).sigma_alpha = sigma_gamma(j);
        end
    end

    % Prey consumption
    for i=1:length(preypop.individuals)
        preypop.individuals(i).a_k = preypop.individuals(i).consumption(resources,1);
    end

    preyind = [preypop.individuals];
    cellarr = {preyind.a_k};
    a_k_matrix = cat(3,cellarr{:});
    % Calculate equilibrium resource abundance
    for j=1:3
        % Calc Rk_eq for the kth resource
        resources(j).Rk_eq = resources(j).eq_abundance(2, a_k_matrix(1,j,:));
    end

    % Pred consumption
    for i=1:length(predpop.individuals)
        predpop.individuals(i).a_k = predpop.individuals(i).consumption([preypop.individuals.alpha], 1);
        predpop.individuals(i).fitness = predpop.individuals(i).calc_fitness_alpha(1);
    end
    
    predpop.attack_rate = predpop.update_attack_rate();
    pred_attack = [predpop.attack_rate];
    for i=1:length(preypop.individuals)
        predatt = pred_attack(:,i);
        preypop.individuals(i).fitness = preypop.individuals(i).calc_fitness_alpha(resources,predatt,1);
    end
    
    sigmapred = [predpop.individuals.alpha];
    fitpred = [predpop.individuals.fitness]/2;
    fit = reshape(fitpred, [length(sigma_gamma), length(gamma)]);
    surf(gamma,sigma_gamma, fit);
    ylabel('\sigma_\gamma');
    xlabel('\gamma');
    zlabel('fitness');

end
