%% Side scan sonar simulator
%% Tips


%% Reset system
clc;
clear;
close all;

c = 1500;
%% [A] Engineering components
% (1) Working signal 
FC = 400*10^3; % center frequency
FR = 4*FC;     % refresh frequency
BW = 20*10^3;  % bandwidth
PW = 5*10^-3;  % pulse width
DR = 250;      % detection range
SP = DR*2/c;   % signal period
save('START.mat','FC','FR','BW','PW','DR','SP');
WorkingSignal;
disp('WorkingSignal.m is finished');

% (2) ULA designer
RElement_Num = 179*5; % Number of one single receiving array elements    [Odd number]
RSubarray_Num = 5;  % Number of receiving arrays            [Odd number]
TElement_Num = 79;  % Number of Transmitting array elements [Odd number]
TSubarray_Pos = 2;  % Position of transmitting array [Align receiving array position]
ULA_InstallationAngle = 20; % Vertical installation angle
ULA_VerticalAngle = 50;     % Vertical beam angle
save('START.mat','RElement_Num','RSubarray_Num','TElement_Num','TSubarray_Pos',...
    'ULA_InstallationAngle','ULA_VerticalAngle','-append');
ULADesigner;
ULA_BeamPattern;
disp('ULADesigner.m is finished');

% (3) Navigation trajectory
load('START.mat','SP');
% ① Speed
NavigationSpeed_knot = 2;                      % (mile/h)
NavigationSpeed_mpers = NavigationSpeed_knot/2; % (m/s)
% ② Displacement
Voyage = 20;
Tracks = 0:NavigationSpeed_mpers*SP:Voyage+NavigationSpeed_mpers*SP; % forward is +
Sway = zeros(1,length(Tracks)); % right is +
Heave= zeros(1,length(Tracks)); % down is +
% Sway = real(exp(1i*(2*pi*(1/7))*(1:length(Tracks))));
% Heave = 2*real(exp(1i*(2*pi*(1/7))*(1:length(Tracks))));
% ③ Attitude
Roll = zeros(1,length(Tracks));  % 
Yaw = zeros(1,length(Tracks));   % right is +
Pitch = zeros(1,length(Tracks)); % down is +
save('START.mat','NavigationSpeed_mpers','Voyage','Tracks','Sway','Heave',...
    'Roll','Yaw','Pitch','-append');
NavigationTrajectory;
disp('NavigationTrajectory.m is finished');

%% [B] Aerial survey simulation component
load('START.mat','Tracks','RSubarray_Num');
TRANSMITTER = 1:length(Tracks)-1; % 1:length(Tracks)-1
RSUBARRAY = 2;  % 1:RSubarray_Num
save('START.mat','RSUBARRAY','TRANSMITTER','-append');

% (1) Seabed model(Built in settings)
% Notes: different experimental purposes are to use different seabed models
SeabedGenerator_Curve;
% SeabedGenerator_KI;
% TMGenerator_KI;
% SeabedGenerator_Cylindrical;
% SeabedGenerator_Lake;
disp('SeabedGenerate.m is finished');

% (2) Navigation scanning
Scoring_criteria = 1*0; % simplification percentage   0.65 is right
save('START.mat','Scoring_criteria','-append');
NavigationScanning;
disp('NavigationScanning.m is finished'); 

% (3) Sonar parameters
SonarParameters;
disp('SonarParameters.m is finished');

% (4) Array signal synthesis
load('START.mat','FC');
FS = 4*FC; % AD sampling frequency
save('START.mat','FS','-append');
ArraySignalSynthesis;
disp('ArraySignalSynthesis.m is finished');

% (5) Analog circuit processing
noise_fl = 370*10^3; % noise frequency
noise_fh = 430*10^3; % noise frequency
Fixed_gain_dB = 0;
save('START.mat','noise_fl','noise_fh','Fixed_gain_dB','-append');
AnalogCircuitProcessing;
disp('MBSSS_AnalogCircuitProcessing.m is finished');
%% [C] Algorithm processing module
% [1] Single-beam side-scan sonar
% (1-1) SBSSS Digital signal processing
SBSSS_DigitalSignalProcessing;
disp('SBSSS_DigitalSignalProcessing.m is finished');
% (1-2) SBSSS imaging
SBSSS_ImageDisplay;
disp('SBSSS_ImageDisplay.m is finished');


% % % % % [2] Multi-beam side-scan sonar
% % % % % (2-1) Synthetic receiving aperture
% % % % MBSSS_SyntheticAperture;
% % % % disp('MBSSS_SyntheticReceivingAperture.m is finished');
% % % % % (2-2) Array manifold
% % % % MBSSS_ArrayManiFold;  
% % % % disp('MBSSS_ArrayManiFold.m is finished');
% % % % % (2-3) Beam forming processing(include digital signal processing)
% % % % MBSSS_Beamforming;
% % % % disp('MBSSS_Beamforming.m is finished');
% % % % % (2-4) MBSSS imaging
% % % % MBSSS_ImageDisplay;
% % % % disp('MBSSS_ImageDisplay.m is finished');












