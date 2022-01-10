function deterministic_main
    t_end = 2000;
    N_vec = zeros(1,t_end);
    P_vec = zeros(1,t_end);
    R_vec = zeros(1,t_end);
    fn_vec = zeros(1,t_end);
    fp_vec = zeros(1,t_end);
    N_vec(1) = 100; % N init
    P_vec(1) = 5; % P init
    fn_vec(1) = 2; % N init
    fp_vec(1) = 2; % P init
    F = 2; % Fecundity rate
    K = 500; % Carrying capacity
    a0N = 3; % Base attack rate prey
    a0P = 0.008; % Base attack rate pred
    bmax = 1; % Pred eff
    g = 0.5; % Conversion
    outfile = fopen('parameters.csv','w');
    fprintf(outfile,'%f\t%f\t%f\t%f\t%f\n', ["a_N", "a_P", "g", "N", "P"]);
    for a0N=1:0.1:4
        for a0P=0:0.001:0.5
            for g=0:0.01:1
                for i=2:t_end
                    R = resource_ab(K,N_vec(i-1),a0N,F);
                    R_vec(i) = R;
                    fitness_pred = fit_pred(g,bmax, N_vec(i-1), a0P);
                    fitness_prey = fit_prey(a0N, R, K, bmax, P_vec(i-1), a0P);
                    N_1 = prey_fun(N_vec(i-1), F, fitness_prey);
                    if N_1>0 N_vec(i) = N_1; else N_vec(i) = 0; end;
                    P_1 = pred_fun(P_vec(i-1), F, fitness_pred);
                    if P_1>0 P_vec(i) = P_1; else P_vec(i) = 0; end;
                end
                if N_1>1 
                    fprintf(outfile,'%f\t%f\t%f\t%f\t%f\n', [a0N, a0P, g, N_1, P_1]);
                end
            end
        end
    end
    fclose(outfile);
    figure;
    plot(1:t_end, N_vec);
    hold on
    plot(1:t_end, P_vec);
    legend(["N", "P"]);
    hold off
    figure;
    plot(2:t_end, fn_vec(2:t_end));
    hold on
    plot(2:t_end, R_vec(2:t_end));
    hold off
end

function N_1 = prey_fun(N, F, fitness)
    N_1 = N * fitness * F/F;
end

function P_1 = pred_fun(P, F, fitness)
    P_1 = P * (fitness/F) * F;
end

function R = resource_ab(K, N, a0N, F)
    R = K/(1+(N*a0N/(K*F)));
end

function fi = fit_prey(a0N, R, K, bmax, P, a0P)
    fi = (R*a0N/K) - (bmax*a0P*P);
end

function fi = fit_pred(g, bmax, N, a0P)
    fi = g*bmax*a0P*N;
end