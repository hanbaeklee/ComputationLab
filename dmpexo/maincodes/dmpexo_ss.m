%% A canonical DMP model with exogenous separation 
%April. 2023
%steady-state
%=========================
%housekeeping
%=========================
clc;
clear variables;
close all;
fnpath = '../functions';
addpath(fnpath);

%=========================
%parameter setup
%=========================
pa.m        = 0.2910;    % Match efficiency
pa.xi1      = 0.6353;    % Matching function parameter CD
pa.xi2      = 0.5250;    % Matching function parameter CES
pa.beta     = 0.9900^(1/3); % Discount factor
pa.sigma    = 1.0000;    % Intertemporal elasticity of substitution
pa.lambda   = 0.0283;    % Exogenous separations
pa.kappa    = 0.0667;    % Vacancy posting cost
pa.b        = 0.6720;    % Unemployment benefit
pa.eta      = pa.xi1;    % Hosios (1990) Condition for CD (can write it for CES later)
% pa.eta      = pa.xi1;   % Hosios (1990) Condition for CD (can write it for CES later)

%%
%=========================
%stead-state equilibrium
%=========================
error       = 10;
weightOld   = 0.9;
A           = 1;
% w           = (1-pa.eta)*pa.b + pa.eta*1;
q           = 0.8;

while error > 1e-10

v           = (q/pa.m)^(1/(-pa.xi1)) / ( 1+(1/pa.lambda)*q*(q/pa.m)^(1/(-pa.xi1)) );
n           = v*q/pa.lambda;

u           = 1-n;
theta       = v/u;
w           = (1-pa.eta)*pa.b + pa.eta*(A + pa.kappa*theta);
c           = A*n - pa.kappa*v + (1-n)*pa.b;

qNew        = pa.kappa/(pa.beta*(A-w+(1-pa.lambda)*pa.kappa/q));
error       = abs(q-qNew);
q           = weightOld*q + (1-weightOld)*qNew;

end

%=========================  
% Report
%=========================  
fprintf(' \n');
fprintf('Unemployment rate: %.10f \n', 1-n);
fprintf('Job finding rate: %.10f \n', q); 

%=========================
%save
%=========================
dir = '../solutions/dmpexo_ss';
save(dir);
