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
Cp = 0.825;   %[kJ/kg*K]
h = 5;  %[W/m^2*K]
T_inf = 298;    %[K]
V_b = 3.42E-5;  %[m^3]
L_b = 65.15E-3; %[m]
A_s = V_b/L_b;  %[m^2]

% Time step
dt = 1; %[hr]

% Tuning parameters;
a_c = 50;  %alpha/cost ratio
V_oc = 3.6; %[volts]
Q_max = 1*3600;  %[coulombs/hr]
Ah = 2.2; %amp-hour throughput

% Limits
P_batt_max = 50; %[kW]
P_dem_max = 2E3; 
G_max = 3E3; %Taken from CE 295 HW3

E_min = 0.3*Ah*V_oc/1000;  %[kWh]
E_max = 0.9*Ah*V_oc/1000;
T_min = 26+273; %[K]
T_max = 36+273; 

% Data
num_nodes = 3; %[3 nodes]
N = 24; %[hr]

t = linspace(0,1,N);    %[hr]

% Demand data subroutine - derived from CE 295 HW4
P_dem = demand_data(P_dem_max);
G = G_max;
Wh_max = N*E_max;

%% Grid State and Preallocate
E_grid = (E_min:5:E_max)';
Wh_grid = (0:5:Wh_max)';

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
    for idx=1:ns
            % Find dominant bounds
            lb = max([-P_dem(k), -P_batt_max, (-E_grid(idx)+E_min)/dt]);
            ub = min([(-E_grid(idx)+E_max)/dt, P_batt_max, G-P_dem(k)]);
        
            % Grid Battery Power between dominant bounds
            P_batt_grid = linspace(lb,ub,200)';
            
            % Calculate log of capacity loss
            Q_k = log(B) + E_a/(R*T_inf) + z.*log(Wh(k+1)/V_oc);
            
            % Cost-per-time-step (vectorized for all P_batt_grid)
            g_k = a_c*Q_k + P_batt_grid + P_dem(k);
        
            % Calculate next Wh (vectorized for all P_batt_grid)
            E_nxt = E_grid(idx) + dt.*P_batt_grid;
        
            % Compute value function at nxt time step (need to interpolate)
            V_nxt = interp1(E_grid,V(:,k+1),E_nxt,'linear');
        
            % Value Function (Principle of Optimality)
            [V(idx, k), ind] = min(g_k + V_nxt);
        
            % Save Optimal Control
            u_star(idx,k) = P_batt_grid(ind);        
    end
end

solveTime = toc;
fprintf(1,'DP Solver Time %2.2f sec \n',solveTime);

%% Simulate Results

% Preallocate
E_sim = zeros(N,1);
P_batt_sim = zeros(N,1);
Q_sim = zeros(N,1);
Q_cap = zeros(N,1);

% Initialize
E_o = E_max;    
E_sim(1) = E_o;
Wh_sim = 0;

% Simulate Grid Dynamics
for k = 1:(N)
    
    % Use optimal battery power, for given Wh (need to interpolate)
    P_batt_sim(k) = interp1(E_grid,u_star(:,k),E_sim(k),'linear');
    
    % Capacity loss
    Q_sim(k) = log(B) + E_a/(R*T_inf) + z*log(Wh_sim/V_oc);
    
    % Capacity
    Q_cap(k) = Q_max*(1-Q_sim(k)/100);
    
    % Time-step battery dynamics
    Wh_sim = Wh_sim + abs(P_batt_sim(k))*dt;
    
    % Energy state evolution
    E_sim(k+1) = E_sim(k) + P_batt_sim(k)*dt;
    
end

fprintf(1,'Total throughput %2.2f kWh\n',Wh_sim);
fprintf(1,'Total Capacity Fade %2.2f C \n',sum(exp(Q_sim/100))*Q_max);
