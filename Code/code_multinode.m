%% CE 295 - Energy Systems and Control
%   Characterizing cumulative battery aging 
%   Yu-Hsin Huang, Yaser Marafee, Raja Selvakumar

clear; close all; clc;
fs = 15;    % Font Size for plots

%% Parameters and Data

% Battery kinetics
B = 30330;  %[-]
E_a = -31500;   %[J/mol]
R = 8.314;  %[J/mol*K]
z = 0.9;  %[-]
R_b = 13E-3;   %[ohms]

% Temperature dynamics   
rho = 1824; %[kg/m^3]
Cp = 0.825;   %[kJ/kg*K]
h = 5E-3;  %[kW/m^2*K]
T_inf = 298;    %[K]
V_b = 3.42E-5;  %[m^3]
L_b = 65.15E-3; %[m]
A_s = V_b/L_b;  %[m^2]

% Time step
dt = 1; %[h]

% Tuning parameters;
a_c = 50;  %alpha/cost ratio
V_oc = 360; %[volts]
Q_max = 1*3600;  %[coulombs/hr]

% Limits
P_batt_max = 50; %[kW]
P_dem_max = 2E3; 
G_max = 3E3; %Taken from CE 295 HW3

E_min = 0.3*Q_max*V_oc/1000;  %[kWh]
E_max = 0.9*Q_max*V_oc/1000;
Wh_max = P_dem_max*dt;

T_min = 26+273; %[K]
T_max = 36+273; 

% Data
num_nodes = 3;
N = 24; %[hr]

t = linspace(0,1,N);    %[hr]
P_dem = demand_data(P_dem_max);
G = G_max/num_nodes;    %power equally distributed to each node 

%% Grid State and Preallocate

% Construct state grids
num_interpl_pts = 100;
E_grid = linspace(E_min, E_max, num_interpl_pts)';
Wh_grid = linspace(0, Wh_max, num_interpl_pts)';

% Grid size
ns = [length(E_grid), length(Wh_grid)];  % No. of states

% Planning horizon (time steps)
N = length(t);

% Preallocate Value Function (rows index state, columns index time)
V = inf*ones(num_nodes, ns(1), ns(2), N+1);

% Preallocate Control (rows index state, columns index time)
u_star = zeros(num_nodes, ns(1), ns(2), N);

%% Solve DP
tic;

% Boundary Condition of Value Function (Principle of Optimality)
V(:,:,N+1) = 0;

% Iterate backward in time
for k = N:-1:1
    % Iterate over all states      
    for ii=1:ns(1)
        for jj = 1:ns(2)
            for nn = 1:num_nodes
            
            if k == N-1
            % Change upper bound on final state
                lb = max([P_dem(k) - G, -P_batt_max, (E_grid(ii)-E_max)/dt]);
                ub = min([(E_grid(ii)-0.5*E_max)/dt, P_batt_max, P_dem(k)]);
            else
            % Find dominant bounds
                lb = max([P_dem(k) - G, -P_batt_max, (E_grid(ii)-E_max)/dt]);
                ub = min([(E_grid(ii)-E_min)/dt, P_batt_max, P_dem(k)]);
            end
            
            % Grid Battery Power between dominant bounds
            P_batt_grid = linspace(lb,ub,200)';
            
            % Calculate capacity loss
            Q_k = log(B) + E_a/(R*T_inf) + z.*log(Wh_grid(jj)); 
            
            %Cost-per-time-step (vectorized for all P_batt_grid)
            g_k = a_c*exp(Q_k) + (P_dem(k) - P_batt_grid);
             
            %Battery Dynamics
            E_nxt = E_grid(ii) - P_batt_grid*dt;
            Wh_nxt = Wh_grid(jj) + abs(P_batt_grid)*dt;
            
            % Compute value function at nxt time step 
            K_grid = linspace(k, k+1, num_interpl_pts);
            V_nxt = interp3(E_grid,Wh_grid, K_grid, permute(V(nn,:,:,:,:), [3 2 1]),...
                E_nxt,Wh_nxt,k+1,'linear');
        
            % Value Function (Principle of Optimality)
            [V(nn,ii,jj, k), ind] = min(g_k + V_nxt);
        
            % Save Optimal Control
            u_star(nn,ii,jj,k) = P_batt_grid(ind);
            end
        end      
    end
end
solveTime = toc;
fprintf(1,'DP Solver Time %2.2f sec \n',solveTime);

%% Simulate Results

% Preallocate
E_sim = zeros(N,1);
Wh_sim = zeros(N,1);
P_batt_sim = zeros(N,1);
Q_sim = zeros(N,1);

% Initialize   
E_sim(1) = E_min;

% Simulate Grid Dynamics
for k =1:N-1
    
    % Use optimal battery power (need to interpolate)
    P_batt_sim(k) = interp2(E_grid,Wh_grid, u_star(:,:,k)', E_sim(k), Wh_sim(k), 'linear');
        
    %  Dynamics
    Wh_sim(k+1) = Wh_sim(k) + abs(P_batt_sim(k))*dt;
    E_sim(k+1) = E_sim(k) - P_batt_sim(k)*dt;
        
    % Capacity fade
    Q_sim(k) = log(B) + E_a/(R*T_inf) + z*log(Wh_sim(k+1));

end
Pgrid =  P_dem - P_batt_sim;
fprintf(1,'Total throughput %2.2f kWh\n', Wh_sim(N));
fprintf(1,'Percentage Capacity Fade %2.2f  \n',(sum(exp(Q_sim))));

%% Plotting
plot(t, P_batt_sim, t ,Pgrid, t, P_dem, t, E_sim, 'linewidth', 2)
xlabel('Time (sec)','FontSize',fs)
ylabel('Power (W)','FontSize',fs)
title('Grid Dynamics w/ SOC = 50%') 
set(gca,'FontSize',fs)
legend('Pbatt','Pgrid','Pdemand','Battery Capacity')

