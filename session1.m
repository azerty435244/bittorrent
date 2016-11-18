%% Loading MOS tables
addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));
clear;
close all;
clc;

load ('UMC65_RVT.mat'); %Data base for analog circuit

%% Initialize everything
designkitName   = 'umc65'; % minimum gate length 65nm
circuitTitle    = 'Analog Design - Session 1';

%Declaration of the circuit components
elementList.nmos = {'Mn1','Mn2','Mn3'};
elementList.pmos = {};
 
spec.VDD        = 1.1;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator       ='spectre';
simulFile       = 0;
simulSkelFile   = 0;
analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice,...
    designkitName, NRVT, PRVT, simulator, simulFile, simulSkelFile);

analog          = cirCheckInChoice(analog, choice);

%% 
disp('=====================================================');
disp('= Session 1: Transistor curves + Matlab environment =');
disp('=====================================================');
% Welcome to the lab sessions of the course of Analog Electronics. This
% matlab file is a starting point for the lab. As you can notice, the
% exercices are guided with some extra information and some questions
% marked with "Q:". Understanding the concepts, discussed during the labs 
% (and theory), must prepare you to pass the exam. Extra exercises are
% given as homework. This is not mandatory, but can be used to train your
% skills further. Homework is marked with "HW:". 

% This lab session consists of three exercises. Most of the code is already 
% provided in the proper section. You only have to implement code between 
% the 'Write your code from here' and the 'Until here' brackets. If you
% store the data in the correct vectors, the provided code will plot for
% you, so use the same vector names!!

% Good luck!!

% TIP: if you don't know how to implement a function in matlab, use the
% 'doc' command to get more information from matlab. For example to know
% more about the 'plot' function, type in the Command Window: doc plot
%%%%%%%%%%%%%%%%%%%%%%%%%%% Circuit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('        Vd            ');
disp('        |             ');
disp('  Vg---Mnx            ');
disp('        |             ');
disp('        Vs            ');

fprintf('\n--- First Exercise: Designing transistor from scratch ---\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Goal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example exercise:
% Play with VD, VS, VOV, L and W and try to make expectations about gm, IDS
% and the different capacitances. For example: increase VOV with a factor 2
% or decrease W with a factor 2. How do the parameters change??
% Look as well to the syntax of the exercise. Make yourself use to the
% functions used

%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VD      = spec.VDD/2;                   % [V] Drain voltage
VS      = 0;                            % [V] Source voltage

Mn1.lg  = 200e-9;                       % [m] Design choice
Mn1.vov = 0.40;                         % [V] Design choice: Choice of inversion
Mn1.vsb = 0;                            % [V] Bulk connected to lowest potential
Mn1.vds = VD-VS;                        % [V] drain source voltage
Mn1.vth = tableValueWref('vth',NRVT,Mn1.lg,0,Mn1.vds,Mn1.vsb); % [V] look up the treshold voltage given these conditions
Mn1.vgs = Mn1.vov + Mn1.vth;            % [V] Gate source voltage: By definition
Mn1.w   = 1000e-9;                      % [m] width of the transistor
Mn1     = mosNfingers(Mn1);             % [] Design choice: cut total width in pieces
Mn1     = mosOpValues(Mn1);             % Calculating operating point values
% such as DC current, small signal parameters such as gds, parasitic
% capacitances, noise figure

%%%%%%%%%%%%%%%%%%%%% Figures of Merit + Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nTransistor parameters:\n');
fprintf('W = %gnm\nL = %gnm\nVOV = %6.2fV\nVDS = %6.2fV\n\n', Mn1.w*1e9, Mn1.lg*1e9, Mn1.vov, Mn1.vds);
fprintf('gm = %gmMho\ngds = %gmMho\ngm/gds = %g\nIDS = %gmA\n',Mn1.gm*1e3,Mn1.gds*1e3,Mn1.gm/Mn1.gds,Mn1.ids*1e3);
fprintf('CGS = %gfF\nCDS = %gfF\nCGD = %gfF\nfT = %gMHz\n', Mn1.cgs*1e15, Mn1.cdd*1e15,Mn1.cgd*1e15,Mn1.ft/1e6);

% Q: What happens with gm, IDS when you double the width while keeping the 
% other parameters constant?
% Q: How do W and VOV scale, if you want to make gm 4 times as large while 
% keeping the same current consumption? To which inversion level are moving
% to?
% Q: How do W and VOV scale, if you want to keep gm constant while reducing
% power consumption with a factor 4
% gm = cte * W/L * VOV = sqrt(cte * W/L * IDS) = (2 * IDS)/VOV
% Q: A) Does gm scales linearly or with the square root of the current
% B) In what situation do you use which equation?
% C) In which inversion level are these equation applicable

%%
fprintf('\n--- Second Exercise: sweeping IDS and gm/IDS in function of VGS ---\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Goal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this exercise, we are going to sweep 'VGS' and look how the current
% 'IDS' and 'gm/IDS' change while keeping the other paramters constant

% Recycle the code from above to sweep IDS and gm/IDS in function of VGS
% Mn2.vgs has to be swept from 0 till spec.VDD
% !!! Store the values of Mn2.vgs in a vector VGS, store the current
% Mn2.ids in vector IDS, and store the gm/IDS Mn3.gm/Mn2.ids in vector gmIDS
% Remark: You don't have to calculate the threshold voltage in this step
% since you don't need vov or vth to calculate the operating point values.
% To calculate the operating point values, you will need 'VGS' instead.

% TIP: You are working here with Mn2, NOT Mn1. Change the index!!

%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Write your code from here





%// Until here

%%%%%%%%%%%%%%%%%%%%% Figures of Merit + Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(221); plot(VGS,IDS*1e3,'x-'); 
title('Current of Mn2 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('I_{DS} [mA]');
grid on;
subplot(222); semilogy(VGS,IDS*1e3,'x-'); 
title('Current of Mn2 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('I_{DS} [mA]');
grid on;
subplot(223); plot(VGS,gmIDS,'x-'); 
title('g_m/I_{DS} of Mn2 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('g_m/I_{DS} [1/V]');
grid on;
subplot(224); semilogy(VGS,gmIDS,'x-'); 
title('g_m/I_{DS} of Mn2 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('g_m/I_{DS} [1/V]');
grid on;

% Q: Does the transistor current seem to scale linearly with VOV?
% Q: In which inversion level is the transistor most power efficient?

%%
fprintf('\n--- Third Exercise: Designing W in function of VGS for constant current ---\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Goal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this exercise we are going to sweep 'VGS'. We will use a new function
% called mosWidth to calculate the width such that the transistor consumes
% the current we want.

% Recycle the code from above to design 'W' in function of 'VGS' for constant
% current of 0.5mA. 'Mn3.vgs' has to be swept from  till 'spec.VDD'
% Store the values of 'Mn3.vgs' in a vector 'VGS', and store the Width
% 'Mn3.w' in vector 'W'
% At the same time save the cut off frequency of the transistor 'Mn3.ft' in
% the vector 'FT'
% For this exercise you will use function 'mosWidth', this function will give
% back the width of a transistor needed for a given gm or current

%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Write your code from here


% Mn3.w   = mosWidth('ids',Mn3.ids,Mn3);  % desiging width with respect to current


%// Until here

%%%%%%%%%%%%%%%%%%%%% Figures of Merit + Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(221); plot(VGS,W,'x-'); 
title('Width of Mn3 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('W [m]');
grid on;
subplot(222); semilogy(VGS,W,'x-'); 
title('Width of Mn3 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('W [m]');
grid on;
subplot(223); plot(VGS,FT,'x-'); 
title('Cutoff frequency of Mn3 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('Cut off frequency f_T [m]');
grid on;
subplot(224); semilogy(VGS,FT,'x-'); 
title('Cutoff frequency of Mn3 vs gate-source voltage');
xlabel('VGS [V]');
ylabel('Cut off frequency f_T [m]');
grid on;

% Q: What is the disadvantage of going in weak inversion?
% Q: In which inversion level is the transistor fastest?
% HW: Plot fT Mn2.ft in function of L. Does the transistor gets faster:
% A) with smaller L or B) with bigger L, why?