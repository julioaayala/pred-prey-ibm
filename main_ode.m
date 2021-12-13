function main_ode(varargin)
     parameters.a_0N = 1;
     parameters.a_0P = 0.3;
     parameters.d = 0.6;
     parameters.g = 0.5;
     parameters.K = 50;
     parameters.F = 2;
     t_end = 200;
     t = 1:t_end;
     y0 = [10,10];
     [t,y] = ode45(@(t,y) pred_prey_ode(t,y,parameters), t, y0);
     
     figure
     plot(t,y(:,1),'-',t,y(:,2),'-');
     hold on;
     legend('Prey','Predator');
     xlabel('Time');
     ylabel('Population');
     hold off;
 end