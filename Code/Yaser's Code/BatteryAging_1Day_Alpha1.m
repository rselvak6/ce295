%% CE 295 - Energy Systems and Control
%   Characterizing cumulative battery aging 
%   Yu-Hsin Huang, Yaser Marafee, Raja Selvakumar

clear; close all; clc;
fs = 15;    % Font Size for plots

%% Parameters and Data

% Battery kinetics
B = 5000;  %[-]
E_a = -31500;   %[J/mol]
R = 8.314;  %[J/mol*K]
z = 0.55;  %[-]
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
days = 1;

% Tuning parameters;
alpha = 50;  %alpha/cost ratio
BatteryCap = 2.2*3.6e4; %Wh (Ah * V)
% Limits

E_min = 0.01*BatteryCap;  %[Wh]
E_max = BatteryCap;
T_min = 26+273; %[K]
T_max = 36+273; 

P_batt_max = 0.9*E_max/dt; %[W]
P_dem_max = 25E3; 
G_max = 50E3;
Wh_max = 24*days*E_max;

% Data

N = 24*days; %[hr]

t = linspace(0,1,N);    %[hr]
P_dem0 = P_dem_max.*[0.08 0.15 0.21 0.29...
    0.33 0.48 0.51 0.55 ...
    0.84 0.94 0.96 0.72...
    0.45 0.35 0.35 0.39...
    0.63 0.79 0.83 0.82...
    0.81 0.76 0.54 0]'; 
P_dem = repmat(P_dem0, days,1);

cost_k0 = 100*[0.04 0.05 0.07 0.09...
    0.11 0.13 0.16 0.25 ...
    0.48 0.81 0.94 0.87...
    0.77 0.75 0.76 0.74...
    0.87 0.81 0.91 0.94...
    0.85 0.81 0.74 0]'; 
cost_k = repmat(cost_k0, days,1);

%Assume each node draw the same maximum amount of power from the grid
G = G_max;

%% Grid State and Preallocate
E_grid = linspace(E_min, E_max, 50)';
Wh_grid = linspace(0, Wh_max, 50)';

% Grid size
ns = [length(E_grid), length(Wh_grid)];  % No. of states

% Planning horizon (time steps)
N = length(t);

% Preallocate Value Function (rows index state, columns index time)
V = inf*ones(ns(1), ns(2), N+1);

% Preallocate Control (rows index state, columns index time)
u_star = zeros(ns(1), ns(2), N);

%% Solve DP
tic;

% Boundary Condition of Value Function (Principle of Optimality)
V(:,:,N+1) = 0;

% Iterate backward in time
for k = N:-1:1
    % Iterate over all states
        
    for idx=1:ns(1)
        for ii = 1:ns(2)
            
            if k == N-1
            % Find dominant bounds
            lb = max([P_dem(k) - G, -P_batt_max, ... 
                (E_grid(idx)-E_max)/dt]);
            ub = min([(E_grid(idx)-0.3*E_max)/dt, ...
                P_batt_max, P_dem(k)]);
            else
            % Find dominant bounds
            lb = max([P_dem(k) - G, -P_batt_max, ... 
                (E_grid(idx)-E_max)/dt]);
            ub = min([(E_grid(idx)-E_min)/dt, ...
                P_batt_max, P_dem(k)]);
        
            end
            
        
            % Grid Battery Power between dominant bounds
            P_batt_grid = linspace(lb,ub,200)';
            
            
            % Calculate capacity loss
            Q_k = log(B) + E_a/(R*T_inf) + z.*log(Wh_grid(ii)*2.2/E_max);
            %Scalling Wh as it was for a 2.2Ah battery
            
            %Cost-per-time-step (vectorized for all P_batt_grid)
            g_k = alpha*E_max*exp(Q_k) + cost_k(k).*(P_dem(k) - P_batt_grid);
             
            %Battery Dynamics
            E_nxt = E_grid(idx) - P_batt_grid*dt;
            Wh_next = Wh_grid(ii) + abs(P_batt_grid*dt);
            
            % Compute value function at nxt time step 
            V_nxt = interp2(E_grid,Wh_grid, V(:,:,k+1)',E_nxt,Wh_next,'linear');
        
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
        
    %  Dynamics
    Wh_sim(k+1) = Wh_sim(k) + abs(P_batt_sim(k))*dt;
    E_sim(k+1) = E_sim(k) - P_batt_sim(k)*dt;
        
    % Capacity fade
    Q_sim(k) = log(B) + E_a/(R*T_inf) + z*log(Wh_sim(k+1)*2.2/E_max);
    
    % Time-step battery dynamicsE_sim(k+1) = E_sim(k) - dt*P_batt_sim(k);
end
Pgrid =  P_dem - P_batt_sim;
time = linspace(1,N,N);
fprintf(1,'Total throughput %2.2f Wh\n', Wh_sim(N));
fprintf(1,'Percentage Capacity Fade %2.2f  \n',(sum(exp(Q_sim))));

%% Plotting
subplot(3,1,1)
hold on
plot(time, P_batt_sim/1000, 'blue', 'linewidth', 2)
plot(time, Pgrid/1000,'r:', 'linewidth', 2)
plot(time, P_dem/1000,'black','linewidth', 1);
hold off
xlabel('Time [hr]','FontSize',10,'FontWeight','bold');
ylabel('Power [kW]','FontSize',10,'FontWeight','bold');
title('$$\alpha$$ = 1','Interpreter','latex','FontSize',20);
legend('P_{batt}','P_{grid}','P_{dem}')
set(gca,'Fontsize',fs)

subplot(3,1,2)
for ij = 1:N
    if ij == 1
       Q_sim_cum(ij) = exp(Q_sim(ij));
    else
       Q_sim_cum(ij) = exp(Q_sim(ij)) + Q_sim_cum(ij-1);
    end 
end
plot(Wh_sim/1000, (100-Q_sim_cum)/100,'blue','linewidth',2)
xlabel('Power throughput [kWh]','FontSize',10, 'FontWeight','bold');
ylabel('Capacity [%]','FontSize',10,'FontWeight','bold');
set(gca,'Fontsize',fs)
Emin = E_min/E_max*ones(length(time),1);
Emax = ones(length(time),1);

subplot(3,1,3)
plot(time, E_sim/E_max, 'blue', 'linewidth',2)
hold on
plot(time,Emin,'r--','linewidth',1)
plot(time,Emax,'r--','linewidth',1);
xlabel('Time [hr]','FontSize',10,'FontWeight','bold');
ylabel('SOC','FontSize',10,'FontWeight','bold');
set(gca,'Fontsize',fs)
hold off