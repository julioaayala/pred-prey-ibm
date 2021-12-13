function deterministic_main
    t_end = 1000;
    N_vec = zeros(1,t_end);
    P_vec = zeros(1,t_end);
    R_vec = zeros(1,t_end);
    fn_vec = zeros(1,t_end);
    fp_vec = zeros(1,t_end);
    N_vec(1) = 365; % N init
    P_vec(1) = 50; % P init
    fn_vec(1) = 2; % N init
    fp_vec(1) = 2; % P init
    F = 2; % Fecundity rate
    K = 1000; % Carrying capacity
    a0N = 2; % Base attack rate prey
    a0P = 0.002; % Base attack rate pred
    bmax = 1; % Pred eff
    g = 0.5; % Conversion

    for i=2:t_end
        R = resource_ab(K,N_vec(i-1),a0N,F);
        R_vec(i) = R;
        fitness_pred = fit_pred(g,bmax, N_vec(i-1), a0P);
        fitness_prey = fit_prey(a0N, R, K, bmax, P_vec(i-1), a0P);
        N_1 = prey_fun(N_vec(i-1), F, fitness_prey);
        if N_1>0 N_vec(i) = N_1; else N_vec(i) = 0; end;
        P_1 = pred_fun(P_vec(i-1), F, fitness_pred);
        if P_1>0 P_vec(i) = P_1; else P_vec(i) = 0; end;
        fn_vec(i) = fitness_prey;
        fp_vec(i) = fitness_pred;
    end
    figure;
    plot(1:t_end, N_vec);
    hold on
    plot(1:t_end, P_vec);
    legend(["N", "P"]);
    hold off
    figure;
    plot(2:t_end, fn_vec(2:t_end));
    hold on
    plot(2:t_end, fp_vec(2:t_end));
    hold off
    disp(N_vec);
    disp(P_vec);
    disp(fn_vec);
    disp(fp_vec);
end

function N_1 = prey_fun(N, F, fitness)
    N_1 = N * fitness * F;
end

function P_1 = pred_fun(P, F, fitness)
    P_1 = P * fitness * F;
end

function R = resource_ab(K, N, a0N, F)
    R = K/(1+(N*(a0N/K)/1));
end

function fi = fit_prey(a0N, R, K, bmax, P, a0P)
    fi = (R*a0N/K) - (bmax*a0P*P);
end

function fi = fit_pred(g, bmax, N, a0P)
    fi = g*bmax*a0P*N;
end