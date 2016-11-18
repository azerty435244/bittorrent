%% Loading MOS tables
addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));
clear;
close all;
clc;

load ('UMC65_RVT.mat');

%% Initialize everything
designkitName   = 'umc65';
circuitTitle    = 'Analog Design - Session 2';

%Declaration of the circuit components
elementList.nmos = {'Mn1','Mn2'};
elementList.pmos = {'Mp2'};
 
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
disp('============================');
disp('= Session 2: Common source =');
disp('============================');

fprintf('\n--- First Exercise: Common source amplifier with resistor ---\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%% Circuit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('       VDD           ');
disp('        |            ');
disp('        R            ');
disp('        |----+-OUT   ');
disp('  IN---Mn1   |       ');
disp('        |    CL      ');
disp('        |    |       ');
disp('       GND  GND      ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec.fGBW       = 100e6;    % [Hz] GBW frequency GainBandwidth
spec.Cl         = 10e-12;   % [F] load capacitance
spec.VDD        = 1.1;      % [V] Power supply voltage (1.1 -> 2.1V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Goal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get at least a gain of 18dB (width of transistor < 1mm, length of transistor <= 200nm)
% Notice that at the end mosCheckSaturation is used to check whether the
% transistor of interest is biased in saturation!

%%%%%%%%%%%%%%%%%%%%%%% Design Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mn1.lg  = 200e-09;               % [m] Design choice (channel length)
Mn1.vov = -0.05;                 % [V] Design choice: Choice of inversion

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Don't touch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mn1.vsb = 0;                            % Bulk connected to lowest potential
Mn1.vds = spec.VDD/2;                   % Maximum swing. Why VDD/2? Biais point on the middle 
% of VDD -> 0 to have the maximum swing
% level: lower vov is more gain, lower current, very big W
Mn1.vth = tableValueWref('vth',NRVT,Mn1.lg,0,Mn1.vds,Mn1.vsb);
Mn1.vgs = Mn1.vov + Mn1.vth;            % By definition
Mn1.gm  = 2*pi*spec.fGBW*spec.Cl;       % gm is fixed by GBW : GBW = gm/(2*pi*Cl)
Mn1.w   = mosWidth('gm',Mn1.gm,Mn1);    % width for given gm. Calculates the width for a given gm
Mn1     = mosNfingers(Mn1);             % Design choice: cut total width in pieces
Mn1     = mosOpValues(Mn1);             % Calculating operating point values
% such as DC current, small signal parameters such as gds, parasitic
% capacitances, noise figure
R       = (spec.VDD - Mn1.vds)/Mn1.ids; % Constrained by the rest of the circuit
Gain    = Mn1.gm/(Mn1.gds+1/R);         
Gain_dB = 20*log10(Gain);              % By definition: Increasing gain is
% increasing gm, decreasing Mn1.gds, increasing R. gm cannot be increased
% since it is set by GBW. You can decrease gds
% by increasing L, or by going into lower inversion since gm/gds is bigger
% in weak inversion. For a given gm, you thus have a smaller gds. Going to
% lower inversion, results in less current through the resistor and
% transistor. As a result R = (VDD- Mn1.vds)/Mn1.ids is bigger. Since R and
% Mn1.gds are increased at the same time, the gain is increased.

% Q: What happens when you increase power supply, when rds>R and when rds<R?
% R: The gain increases more if rds > R
% Q: R>rds when a) in weak inverion with small L, b) strong inversion, big L?
% R: a) in weak inversion(Vov<0) with small L
% Q: What is best? a) rds>R, b) rds = R, c) rds<R
% R: rds>R for maximum gain

% When L is set to L = 100nm instead of Lmin and vov = 0.1V, relative good
% values for the gain can be found. A possibility to increase the gain is
% to increase the voltage of the power supply (eg. 2.1V) since then the 
% voltage over the load resistor is higher for the same current, making it 
% bigger... The drawback is off course a bigger power consumption. Another
% thing we can do is to take a PMOS as a load. A PMOS has a flat
% charactertic around the operating point, while it does not force us to
% use a big power supply.

% Sizing a transistor means looking to the specs, and trying to meet them.
% Therefore trade-offs need to be made. An important thing to remember is
% that you know that when you change operating point and sizes of
% transistors, that the noise behaviour, gain, poles,.. change. And you need
% to get a feeling in what direction they change, eg: increase or decrease 
% of amount of noise, lower frequency of a pole at a certain node...

%%%%%%%%%%%%%%%%%%%%%%% Figures of Merit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nSolution of exercise 1:\n');
pole    = -(Mn1.gds+1/R)/spec.Cl;
freq    = logspace(1,10,1e3);
s       = tf('s');
TF1     = Gain/(1-s/pole);
% figure;
% bode(TF1,freq); grid on;
% title('Voltage gain (common source with resistor load)');

fprintf('R       = %.2f Ohm \n', R);
fprintf('rds     = %.2f Ohm \n', 1/Mn1.gds);
fprintf('Gain    = %.2f \n', Gain);
fprintf('Gain_dB = %.2f dB \n', Gain_dB);
fprintf('Vds     = %.2f V \n', Mn1.vds);
fprintf('Vdsat   = %.2f V \n', Mn1.vdsat);
fprintf('width   = %.2f mm \n', Mn1.w*1000);

fprintf('\nGain(in dB) requirement > 18 dB:\n');
if (Gain_dB > 18)
	fprintf('your design: Success!\n') % Check Gain requirement
else
    fprintf('your design: Fail...\n')
end

fprintf('\nTransisors in saturation:\n');
if mosCheckSaturation(Mn1)
	fprintf('Mn1: Success!\n') % Check if transistor is in saturation
end

fprintf('\nWidth requirement < 1 mm:\n');
if (Mn1.w < 0.001)
	fprintf('Mn1: Success!\n') % Check width requirement
else
    fprintf('Mn1: Fail...\n') 
end

%% 
fprintf('\n--- Second Exercise: Common source amplifier with PMOS ---\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%% Circuit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('       VDD           ');
disp('        |            ');
disp('       Mp2           ');
disp('        |----+-OUT   ');
disp('  IN---Mn2   |       ');
disp('        |    CL      ');
disp('        |    |       ');
disp('       GND  GND      ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec.fGBW       = 100e6;    % [Hz] GBW frequency
spec.Cl         = 10e-12;   % [F] load capacitance
spec.VDD        = 1.1;      % [V] Power supply voltage 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Goal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get at least a gain of 26dB (width of transistor < 1mm, length of transistor <= 200nm)

%%%%%%%%%%%%%%%%%%%%%%% Design Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mn2.lg  = 100e-09;                       % Design choice: flatness
Mn2.vov = 0.1;                       % Design choice: inversion level
Mn2.vds = spec.VDD/2;                       % Design choice: maximum swing
Mp2.lg  = 100e-09;                       % Design choice: flatness
Mp2.vov = -0.1;                       % Design choice: inversion level

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Don't touch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Mn2 instead of reusing Mn1 because nFingers fails when reusing Mn1
Mn2.vsb = 0;                            % For NMOS, bulk is connected to 
% lowest potential
% level: lower vov is more gain, lower current, very big W
Mn2.vth = tableValueWref('vth',NRVT,Mn2.lg,0,Mn2.vds,Mn2.vsb);
Mn2.vgs = Mn2.vov + Mn2.vth;            % By definition
Mn2.gm  = 2*pi*spec.fGBW*spec.Cl;       % gm is fixed by GBW
Mn2.w   = mosWidth('gm',Mn2.gm,Mn2);    % width for given gm
Mn2     = mosNfingers(Mn2);             % Design choice
Mn2     = mosOpValues(Mn2);             % Calculating operating point values

Mp2.vsb = 0;                            % For PMOS, bulk is connected to 
% highest potential
Mp2.vds = Mn2.vds - spec.VDD;           % Forced by choice of vds of NMOS
Mp2.vth = tableValueWref('vth',PRVT,Mp2.lg,0,Mp2.vds,Mp2.vsb);
Mp2.vgs = Mp2.vov + Mp2.vth;            % By definition
Mp2.ids = Mn2.ids;                      % Forced by NMOS
Mp2.w   = mosWidth('ids',Mp2.ids,Mp2);
Mp2     = mosNfingers(Mp2);             % Design choice
Mp2     = mosOpValues(Mp2);             % Calculating operating point values

Gain    = Mn2.gm/(Mn2.gds+Mp2.gds);
Gain_dB = 20*log10(Gain);
%%%%%%%%%%%%%%%%%%%%%%% Figures of Merit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nSolution of exercise 2:\n');
fprintf('rds_n = %6.2f Ohm\n', 1/Mn2.gds);
fprintf('rds_p = %6.2f Ohm\n', 1/Mp2.gds);
fprintf('Gain = %.2f\n',Gain);
fprintf('Gain_dB = %.2f dB\n',Gain_dB);
fprintf('width of Mp2   = %.2f mm \n', Mp2.w*1000);
fprintf('width of Mn2   = %.2f mm \n', Mn2.w*1000);

% At minimum L, the gain is low due to the high gds values of the
% transistors. Increasing length is a way to increase the gain.
% Q: Get the highest gain possible...(ignore the limitation on length but width < 1mm)
% R: 
% Q: Would you put Mp2 in weak inversion?
% R: 

pole    = -(Mn2.gds+Mp2.gds)/spec.Cl;
freq    = logspace(1,10,1e3);
s       = tf('s');
TF2     = Gain/(1-s/pole);
% figure;
% bode(TF2,freq); grid on;
% title('Voltage gain (Common source with PMOS load)');

fprintf('\nGain(in dB) requirement > 26 dB:\n');
if (Gain_dB > 26)
	fprintf('your design: Success!\n') % Check Gain requirement
else
    fprintf('your design: Fail...\n') 
end

fprintf('\nTransisors in saturation:\n');
if mosCheckSaturation(Mp2)
	fprintf('Mp2: Success!\n') % Check if transistor is in saturation
else
    fprintf('Mp2: Fail...\n') 
end
if mosCheckSaturation(Mn2)
	fprintf('Mn2: Success!\n') 
else
    fprintf('Mn2: Fail...\n') 
end

fprintf('\nWidth requirement < 1 mm:\n');
if (Mp2.w < 0.001)
	fprintf('Mp2: Success!\n') % Check width requirement
else
    fprintf('Mp2: Fail...\n') 
end
if (Mn2.w < 0.001)
	fprintf('Mn2: Success!\n') 
else
    fprintf('Mn2: Fail...\n') 
end