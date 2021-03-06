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
z = 0.552;  %[-]
R_b = 13E-3;   %[ohms]

% Temperature dynamics   
rho = 1824; %[kg/m^3]
Cp = 825;   %[J/kg*K]
h = 5;  %[W/m^2*K]
T_inf = 298;    %[K]
V_b = 3.42E-5;  %[m^3]
L_b = 65.15E-3; %[m]
A_s = V_b/L_b;  %[m^2]

% Time step
dt = 1; %[hr]
N = 24; %total simulation time steps

% Tuning parameters
alpha = 50;  %alpha
num_batt = 100;

%Battery values
V_oc = 3.6; %[V]
Ah = 2.2; %[Ah] throughput

% Limits
P_batt_max = Ah*V_oc*num_batt; %[W]
P_dem_max = 1000; 
G_max = 1500; %Taken from CE 295 HW3
E_min = 0.1*Ah*V_oc;  %[Wh]
E_max = 0.9*Ah*V_oc;
Wh_max = N*E_max;
G = G_max;

T_min = 26+273; %[K]
T_max = 36+273; 

% Data
num_nodes = 3; %[3 nodes]
t = linspace(0,1,N);    %[hr]

%% Demand data subroutine - derived from CE 295 HW4
P_dem = demand_data(P_dem_max);
c_k = demand_data(0.12); %currently modelled as proportional to P_dem
a_c = alpha/c_k;

%% Grid State and Preallocate
E_grid = (E_min:0.5:E_max)';
Wh_grid = (0:5:Wh_max)';
a_c = a_c';

% Grid size
ns = [length(E_grid), length(Wh_grid)];  % No. of states

% Preallocate Value Function (rows index state, columns index time)
V = inf*ones(ns(1), ns(2), N+1);

% Preallocate Control (rows index state, columns index time)
u_star = zeros(ns(1), ns(2), N);

%% Solve DP
tic;

% Boundary Condition of Value Function (Principle of Optimality)
V(:,N+1) = 0;

% Iterate backward in time
for k = N:-1:1
    % Iterate over all states
    for idx=1:ns(1)
        for ii=1:ns(2)
            % Find dominant bounds
            lb = max([P_dem(k)-G, -P_batt_max, (E_grid(idx)-E_max)/dt]);
            ub = min([(E_grid(idx)-E_min)/dt, P_batt_max, P_dem(k)]);
        
            % Grid Battery Power between dominant bounds
            P_batt_grid = linspace(lb,ub,200)';
            
            % Calculate log of capacity loss
            Q_k = log(B) + E_a/(R*T_inf) + z.*log(Wh_grid(ii)/E_grid(idx));
            
            % Cost-per-time-step (vectorized for all P_batt_grid)
            g_k = a_c(k)*E_max*exp(Q_k) - P_batt_grid + P_dem(k);
        
            % Calculate next Wh (vectorized for all P_batt_grid)
            E_nxt = E_grid(idx) - P_batt_grid.*dt;
            Wh_nxt = Wh_grid(ii) + abs(P_batt_grid)*dt;
        
            % Compute value function at nxt time step 
            V_nxt = interp2(E_grid,Wh_grid,V(:,:,k+1)',E_nxt,Wh_nxt,'linear');
        
            % Value Function (Principle of Optimality)
            [V(idx,ii, k), ind] = min(g_k + V_nxt);
        
            % Save Optimal Control
            u_star(idx,ii,k) = P_batt_grid(ind);   
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
E_sim(1) = 0.3*E_max;

% Simulate Grid Dynamics
for k =1:N-1
    
    % Use optimal battery power (need to interpolate)
    P_batt_sim(k) = interp2(E_grid,Wh_grid, u_star(:,:,k)', E_sim(k), Wh_sim(k), 'linear');
        
    % Dynamics
    Wh_sim(k+1) = Wh_sim(k) + abs(P_batt_sim(k))*dt;
    E_sim(k+1) = E_sim(k) - P_batt_sim(k)*dt;
        
    % Capacity fade
    Q_sim(k) = log(B) + E_a/(R*T_inf) + z*log(Wh_sim(k)/E_sim(k));
end
P_grid =  P_dem - P_batt_sim;
fprintf(1,'Total throughput %2.2f Wh\n', Wh_sim(N));
fprintf(1,'Percentage Capacity Fade %2.2f%%\n',(sum(exp(Q_sim))));
