function fitness_landscape_test(varargin)
    alpha = 1:0.05:3;
    sigma_alpha = 0.001:0.05:2;
    resources = [Resource(1,1,10000), Resource(2,1,0.1), Resource(3,1,10000)];
    inds = Individual.empty;
    
    
    for j=1:length(sigma_alpha)
        for i=1:length(alpha)
            inds(length(inds)+1) = Individual(alpha(i),1,2,1,sigma_alpha(j),1);
        end
    end
    
    for i=1:length(inds)
        inds(i).a_k = inds(i).consumption(resources,1);
    end
    
    cellarr = {inds.a_k};
    a_k_matrix = cat(3,cellarr{:});
    for i=1:length(resources)
        resources(i).Rk_eq = resources(i).eq_abundance(3, a_k_matrix(:,i,:));
    end
    
    for i=1:length(inds)
        inds(i).fitness = inds(i).calc_fitness_alpha(resources);
    end

    fit = [inds.fitness];
    fit2 = reshape(fit, [length(alpha), length(sigma_alpha)]);
    surf(sigma_alpha, alpha, fit2);
end