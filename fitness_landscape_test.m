function fitness_landscape_test(varargin)
    alpha = 1:0.05:3;
    sigma_alpha = 0.25:0.05:0.50;
    resources = [Resource(1,1,500), Resource(2,1,500), Resource(3,1,500)];
    inds = Prey.empty;
    loci.alpha = 8;
    loci.beta = 8;
    loci.dis = 8;
    loci.pref = 8;
    
    for j=1:length(sigma_alpha)
        for i=1:length(alpha)
            if ismember([2.0], alpha)
                for k=1:1
                    inds(length(inds)+1) = Prey(alpha(i),1,2,1,sigma_alpha(j),1, 1,1,loci,0);
                end
            else
                inds(length(inds)+1) = Prey(alpha(i),1,2,1,sigma_alpha(j),1, 1,1,loci,0);
            end
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
        inds(i).fitness = inds(i).calc_fitness_alpha(resources,0,0);
    end
    
    fit = zeros(1,1);
    for j=1:length(sigma_alpha)
        for i=1:length(alpha)
            fit(i,j) = mean([inds([inds.alpha]==alpha(i) & [inds.sigma_alpha]==sigma_alpha(j)).fitness]);
        end
    end
    
    fit = fit/2;
    %fit2 = reshape(fit, [length(alpha), length(sigma_alpha)]);
    surf(sigma_alpha, alpha, fit);
    xlabel('\sigma_\alpha');
    ylabel('\alpha');
    zlabel('fitness');
end