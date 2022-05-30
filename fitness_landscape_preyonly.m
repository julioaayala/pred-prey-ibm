function fitness_landscape_preyonly(varargin)
    alpha = 1:0.0625:3;
    sigma_alpha = 0.35:0.05:0.55;
    resources = [Resource(1,1,100), Resource(2,1,100), Resource(3,1,100)];
    loci_prey.alpha = 16;
    loci_prey.beta = 8;
    loci_prey.dis = 8;
    loci_prey.pref = 8;

    preypop = Population("prey", 2, 1, 2, 1, 0.45, 1, 1, 0, length(alpha)*length(sigma_alpha), 1e-5, loci_prey);

    for i=1:length(alpha)
      for j=1:length(sigma_alpha)  
        preypop.individuals(((i-1)*length(sigma_alpha))+j).alpha = alpha(i);
        preypop.individuals(((i-1)*length(sigma_alpha))+j).sigma_alpha = sigma_alpha(j);
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

    for i=1:length(preypop.individuals)
        preypop.individuals(i).fitness = preypop.individuals(i).calc_fitness_alpha(resources,0,1);
    end
    
    alphaprey = [preypop.individuals.alpha];
    sigmaprey = [preypop.individuals.sigma_alpha];
    fitprey = [preypop.individuals.fitness]/2;
    fit = reshape(fitprey, [length(sigma_alpha), length(alpha)]);

    surf(alpha,sigma_alpha, fit);
    ylabel('\sigma_\alpha');
    xlabel('\alpha');
    zlabel('fitness');

end
