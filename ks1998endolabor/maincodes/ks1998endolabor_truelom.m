%% An RBC model with heterogeneous firms and irreversible investments
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compute the implied true low of motion
%=========================    
%=========================
% housekeeping
%=========================
clc;
clear variables;
close all;
fnpath = '../functions';
addpath(fnpath);

%=========================
% load the solution
%=========================
load '../solutions/ks1998endolabor_bc.mat';

% prepare the variables of interest
tw      = (1-palpha)*vgridA(tsimpath)'.*(tK(1:end-1)./tsupplyL).^(palpha);
tw      = log(tw);
tK      = log(tK);

%%

R2all = zeros(2,pnumgridA);
MSE2all = zeros(2,pnumgridA);
b1 = zeros(27,pnumgridA);
b2 = zeros(27,pnumgridA);

modelchoice = [2, 3, 4, 5, 6,...
               9, 12, 15, 18, 21, 24, 27];
mtable = zeros(13,7);
mtable2 = zeros(13,7);
for modelselec = 1:length(modelchoice)
for iA = 1:pnumgridA
    
    location = find(tsimpath(1:end-1)==iA);
    Llocation = location-1;
    LLlocation = location-2;
    LLLlocation = location-3;
    LLLLlocation = location-4;
    LLLLLlocation = location-5;
    LLLLLLlocation = location-6;
    LLLLLLLlocation = location-7;
    
    Llocation(Llocation<=1)=1;
    LLlocation(LLlocation<=1)=1;
    LLLlocation(LLLlocation<=1)=1;
    LLLLlocation(LLLLlocation<=1)=1;
    LLLLLlocation(LLLLLlocation<=1)=1;
    LLLLLLlocation(LLLLLLlocation<=1)=1;
    LLLLLLLlocation(LLLLLLLlocation<=1)=1;

    y1 = (tK(location+1));
    y2 = (tw(location));

    x = [ones(length(tK(location)),1) ...
        tK(location) tK(location).^2 tK(location).^3 tK(location).^4 tK(location).^5 ...
        tK(Llocation) tK(Llocation).^2 tK(Llocation).^3 ...
        tK(LLlocation) tK(LLlocation).^2 tK(LLlocation).^3 ...
        tK(LLLlocation) tK(LLLlocation).^2 tK(LLLlocation).^3 ...
        tK(LLLLlocation) tK(LLLLlocation).^2 tK(LLLLlocation).^3 ...
        tK(LLLLLlocation) tK(LLLLLlocation).^2 tK(LLLLLlocation).^3 ...
        tK(LLLLLLlocation) tK(LLLLLLlocation).^2 tK(LLLLLLlocation).^3 ...
        tK(LLLLLLLlocation) tK(LLLLLLLlocation).^2 tK(LLLLLLLlocation).^3 ...
        ];
    x = x(:,1:modelchoice(modelselec));

    [b1temp,~,~,~,R1] = regress(y1,x);
    [b2temp,~,~,~,R2] = regress(y2,x);

    R2all(1,iA) = R1(1);
    MSE2all(1,iA) = R1(4);
    
    R2all(2,iA) = R2(1);
    MSE2all(2,iA) = R2(4);        
    
    if modelselec == length(modelchoice)
    b1(:,iA) = b1temp;
    b2(:,iA) = b2temp;
    end

end

mtable(modelselec+1,4:(3+2*pnumgridA)) = [R2all(1,:),R2all(2,:)];
mtable2(modelselec+1,4:(3+2*pnumgridA)) = [MSE2all(1,:),MSE2all(2,:)];

end

display(R2all);

%% 
% save the coefficients
% save("../solutions/ks1998endolabor_bc_truelomcoeff.mat","b1","b2");

%%
% fill-in the table
mtable(2:end,2:3) = [0,1; 0,2; 0,3; 0,4; 0,5;...
                     1,3; 2,3; 3,3; 4,3; 5,3; 6,3; 7,3 ...
                     ];
mtable2(2:end,2:3) = [0,1; 0,2; 0,3; 0,4; 0,5;...
                     1,3; 2,3; 3,3; 4,3; 5,3; 6,3; 7,3 ...
                     ];
% csvwrite('../solutions/ks1998endolabor_bc_truelomcoeff.csv',mtable,0,0)
