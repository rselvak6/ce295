%% CE 295 - Energy Systems and Control
%   Characterizing cumulative battery aging 
%   Yu-Hsin Huang, Yaser Marafee, Raja Selvakumar

clear; close all;
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
dt = 1; %[h]

% Tuning parameters;
a_c = 0.5;  %alpha/cost ratio
V_oc = 330; %[volts]
Q_max = 1*3600;  %[coulombs]

% Limits
P_dem_max = 30; %[kW]

Wh_min = 1.5*V_oc;  %[J]
Wh_max = 2.5*V_oc;
T_min = 26+273; %[K]
T_max = 36+273; 

% Data
num_nodes = 3; %[3 nodes]
N = 24; %[hr]

t = linspace(0,1,N);    %[hr]
P_dem = P_dem_max*rand(N,num_nodes);   %uniform random distribution, time x states


% Plot Engine efficiency Curve
P_eng = linspace(0,P_eng_max,length(P_dem));
figure(2); clf;
plot(P_eng/1e3, eta_eng(P_eng)) % plot efficiency versus engine power, for total range of engine powers
title('Engine efficiency curve');
xlabel('Engine power [W]');
ylabel('$$\eta$$','Interpreter','latex');

%% Grid State and Preallocate
Wh_grid = (Wh_min:0.005:Wh_max)';

% Grid size
ns = length(Wh_grid);  % No. of states

% Planning horizon (time steps)
N = length(t);

% Preallocate Value Function (rows index state, columns index time)
V = inf*ones(ns,N+1);

% Preallocate Control (rows index state, columns index time)
u_star = zeros(ns,N);

%% Solve DP
tic;

% Boundary Condition of Value Function (Principle of Optimality)
V(:,N+1) = 0;

% Iterate backward in time
for k = N:-1:1
    % Iterate over all states
    for idx=1:ns
        % Iterate over all nodes
        for i=1:num_nodes
            % Find dominant bounds
            lb = max([-P_dem(k,i), -P_batt_max, ...
                (Wh_grid(idx)-Wh_max)/dt]);
            ub = min([(-Wh_grid(idx)+Wh_max)/dt,...
                P_batt_max, G(k)/num_nodes+P_dem(k,i)]);
        
            % Grid Battery Power between dominant bounds
            P_batt_grid = linspace(lb,ub,200)';
            
            % Calculate capacity loss
            Q_k = ln(B) - E_a/(R*T_inf) + z.*ln(Wh_grid(idx)/V_oc); 
            
            % Cost-per-time-step (vectorized for all P_batt_grid)
            g_k = a_c*Q_k + P_batt_grid + P_dem;
        
            % Calculate next Wh (vectorized for all P_batt_grid)
            Wh_nxt = Wh_grid(idx)+ dt.*P_batt_grid;
        
        % Compute value function at nxt time step (need to interpolate)
        V_nxt = interp1(Wh_grid,V(:,k+1),SOC_nxt,'linear');
        
        % Value Function (Principle of Optimality)
        [V(idx, k), ind] = min(g_k + V_nxt);
        
        % Save Optimal Control
        u_star(idx,k) = P_batt_grid(ind);
        end
    end
end
%%

    % Iterate over SOC
    for idx = 1:ns
        
        % Find dominant bounds
        lb = max([-P_batt_max, (Wh_grid(idx)-Wh_max)/dt]);
        ub = min([(-Wh_grid(idx)+Wh_max)/dt,...
            P_batt_max, G(k)/num_nodes+P_dem(k)]);
        
        % Grid Battery Power between dominant bounds
        P_batt_grid = linspace(lb,ub,200)';
        
        % Cost-per-time-step (vectorized for all P_batt_grid)
        g_k = alph*dt.*P_eng./eta_eng(P_eng);
        
        % Calculate next SOC (vectorized for all P_batt_grid)
        SOC_nxt = Wh_grid(idx)+ dt/(Qcap*V_oc).*P_batt_grid;
        
        % Compute value function at nxt time step (need to interpolate)
        V_nxt = interp1(Wh_grid,V(:,k+1),SOC_nxt,'linear');
        
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
SOC_sim = zeros(N,1);
P_batt_sim = zeros(N,1);
P_eng_sim = zeros(N,1);
J_sim = zeros(N,1);

% Initialize
SOC_0 = 0.3;    
SOC_sim(1) = SOC_0;

% Simulate PHEV Dynamics
for k = 1:(N-1)
    
    % Use optimal battery power, for given SOC (need to interpolate)
    P_batt_sim(k) = interp1(Wh_grid,u_star(:,k),SOC_sim(k),'linear');
    
    % Compute engine power
    P_eng_sim(k) = P_dem(k)-P_batt_sim(k);
    
    % Fuel Consumption
    J_sim(k) = alph*P_eng_sim(k)./eta_eng(P_eng_sim(k));
    
    % Time-step SOC dynamics
    SOC_sim(k+1) = SOC_sim(k) + dt/(Qcap*V_oc).*P_batt_sim(k);
    
end

fprintf(1,'Final SOC %2.4f \n',SOC_sim(N));
fprintf(1,'Total Fuel Consumption %2.2f kg \n',sum(J_sim)/1e3);

%% Plot Results
figure(3); clf;

subplot(2,2,1);
% UDDS speed versus time 
plot(t, v_dc);
title('UDDS speed vs. time');
xlabel('Time [s]');
ylabel('Cycle speed [m/s]');
set(gca,'FontSize',fs)

subplot(2,2,2);
% SOC versus time
plot(t, SOC_sim);
title('SOC vs. time');
xlabel('Time [s]');
ylabel('State of charge');
set(gca,'FontSize',fs)

subplot(2,2,3);
% Accumulated fuel consumption [g] versus time
J = zeros(length(J_sim),1);
for i=1:length(J_sim)
    for j=1:i
        J(i) = J(i) + J_sim(j);
    end
end
plot(t, J);
title('Accumulated fuel consumption vs. time');
xlabel('Time [s]');
ylabel('Fuel consumption [g]');
set(gca,'FontSize',fs)

subplot(2,2,4);
% Battery and engine power [kW] versus time
plot(t, P_batt_sim/1e3,t, P_eng_sim/1e3);
title('Engine and battery power vs. time');
xlabel('Time [s]');
ylabel('Power [kW]');
legend('Battery power','Engine power','Location','northwest');
set(gca,'FontSize',fs)
