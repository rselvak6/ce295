%% CE 295 - Energy Systems and Control
%   HW 5 : Optimal Energy Management of PHEV via Dynamic Programming
%   Oski Bear, SID 18681868
%   Prof. Moura
%   Due Apr 25, 2014

%   ENGINE EFFICIENCY CURVE

function out = eta_eng(P_eng)

% polynomial coefficients
p1 =   5.128e-08;
p2 =  -5.927e-06;
p3 =   0.0002652;
p4 =    -0.00583;
p5 =      0.0672;
p6 =   2.622e-05;

% Convert from W to kW
x = P_eng/1e3;

% Compute efficiency curve
out = p1*x.^5 + p2*x.^4 + p3*x.^3 + p4*x.^2 + p5*x + p6;