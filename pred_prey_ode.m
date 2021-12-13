 %Predator-prey function
 function dydt = pred_prey_ode(t,y, parameters)
    dydt = [0;0];
    a_0N = parameters.a_0N;
    a_0P = parameters.a_0P;
    d = parameters.d;  
    g = parameters.g;
    K = parameters.K;
    F = parameters.F;
    dydt(1) = a_0N*y(1)*1/(1 + (y(1)*a_0N/(K*F))) - a_0P*y(1)*y(2);    %Prey
    dydt(2) = g*a_0P*y(1)*y(2) - d*y(2);  %Predator

    % Orig
    %dydt(1) = a_0N*y(1)*(K-y(1))/K - a_0P*y(1)*y(2);    %Prey
    %dydt(2) = g*a_0P*y(1)*y(2) - d*y(2);  %Predator
 end