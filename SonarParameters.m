%% Sonar parameters
%% Tips


%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','NavigationSpeed_mpers','FC','BW','PW','DR','TRANSMITTER','RSUBARRAY');
load('Seabed.mat','Seabed_X','Seabed_Y','Seabed_Zm','Seabed_Z','Seabed_delt_xm','Seabed_delt_ym');
load('ULADesigner.mat','ULA_SubarrayBeamAngle','ULA_Elements');
load('ULA_BeamPattern.mat','Angle','DI_SUBARRAY','DI_ARRAYr');
%% Calculation formula
% Active sonar equation
% SL - 2TL + TS - (NL - DI) - RL = DT
% (1) Echo level √
% EL = SL - 2TL + TS + DI + calibration
% (2) Echo margin
% EM = SL - 2TL + TS - [ NL - ( AG + DI ) ] - RL - DT + calibration
%% Calculation preparation
% (1) Propagation loss coefficient α (dB/km) = 1.093613298*(dB/kyd)
alpha=((0.1*(FC/1000)^2)/(1+(FC/1000)^2)+(40*(FC/1000)^2)/(4100+(FC/1000)^2)...
    +2.75*(10^-4)*(FC/1000)^2+0.003)*1; 
% (2) Backscattering coefficient 【Tips】 Manual adjustment
Backscattering_Flat = -1.262;
Backscattering_Raised = -1.262;

%% Calculation process
% Sound source level (dB)
SL = 220;

for transmitter = TRANSMITTER
% transmitter = 7;
    filename = ['NavigationScanning_',num2str(transmitter),'.mat'];
    load(filename,'TD','RD','IncidenceAngle','DOA');
    
    TL = cell(1,max(RSUBARRAY));
    PT = cell(1,max(RSUBARRAY));
    TS = cell(1,max(RSUBARRAY));
    EL = cell(1,max(RSUBARRAY));
    EM = cell(1,max(RSUBARRAY));
    for Rsubarray = RSUBARRAY
% Rsubarray = 2;
        % (1) Transmission loss
        PropagationDistance_Transmitter = TD{1,Rsubarray};
        PropagationDistance_Receiver = RD{1,Rsubarray};        
        TL_Transmitter = 20*log10(1.093613298*PropagationDistance_Transmitter)+alpha*((1.093613298*PropagationDistance_Transmitter)/1000);   
        TL_Receiver = 20*log10(1.093613298*PropagationDistance_Receiver)+alpha*((1.093613298*PropagationDistance_Receiver)/1000);
        TL{1,Rsubarray} = TL_Transmitter+TL_Receiver;
        
        % (2) Propagation time
        PT{1,Rsubarray} = (PropagationDistance_Transmitter+PropagationDistance_Receiver)/c;
        
        % (3) Target strength
        incidenceAngle = IncidenceAngle{1,Rsubarray};
        % ① Flat seabed
        TargetStrength = 1*((Backscattering_Flat+2*log10((cosd(incidenceAngle)).^2))+10*log10(Seabed_delt_xm*Seabed_delt_ym));     
        % ② Raised seafloor
        TargetStrength(Seabed_Z ~= Seabed_Zm) = 1*(Backscattering_Raised+2*log10((cosd(incidenceAngle(Seabed_Z ~= Seabed_Zm))).^2))+10*log10(Seabed_delt_xm*Seabed_delt_ym); 
        
        TS{1,Rsubarray} = TargetStrength; 
        % (4) Noise level
        S = 3; % Sea state level
        NL = (10*log10((FC/1000)^(-1.7))+6*S+55)+10*log10(BW/1000);          
        
        % (5) Array gain
        AG = 10*log10(sum(ULA_Elements(1:end-1)))+10*log10((2*ULA_SubarrayBeamAngle(1,1))/(c/FC))+3;        
        
        % (6) Directivity index
        horizontalAngle = DOA{1,Rsubarray};
        DI = zeros(size(horizontalAngle));                
% % % %         for yy = 1:size(horizontalAngle,1)
% % % %             for xx = 1:size(horizontalAngle,2)
% % % %                 yy
% % % %                 delt_angle = abs(Angle-horizontalAngle(yy,xx));                
% % % %                 DI(yy,xx) = DI_SUBARRAY(delt_angle == min(delt_angle));
% % % %             end
% % % %         end

        % (7) Detection threshold
        d = 15;
        DT = 10*log10(d/(2*PW))-5*log10(3.981);    
        
        % (8) Reverberation level
        RL = 0;        
        
        % Echo level
        EL{1,Rsubarray} = SL-TL{1,Rsubarray}+TS{1,Rsubarray}+DI-28.923957733255680;
        % Echo margin
        EM{1,Rsubarray} = SL-TL{1,Rsubarray}+TS{1,Rsubarray}-(NL-AG)+DI-RL-DT-28.923957733255680;
        EM_calibration = EM{1,Rsubarray};              
        
        % (6) Display
% % % %         figure(1);        
% % % %         Display = 10.^((EL{1,Rsubarray}-199.4)/20);
% % % %         mesh(Seabed_X(1,:),Seabed_Y(:,1),Display);
% % % %         shading interp;
% % % %         xlabel('Transverse direction（m）','FontSize',15); 
% % % %         ylabel('Navigation direction（m）','FontSize',15);
% % % %         zlabel('Incidence angle（°）','FontSize',15);
% % % % %         set(gca,'ZDir','reverse','FontSize',15);
% % % % %         zlim([15,30]);  
        
        disp(Rsubarray);        
    end
    
    % Save data    
    filename = ['SonarParameters_',num2str(transmitter),'.mat'];
    save(filename,'PT','SL','TL','TS','NL','AG','DI','RL','DT','EL');
    disp(filename);
    sound(sin(2*pi*20*(1:4000)/100));         
end
    
%% Save data
message = [num2str(transmitter),' groups of sonar parameters have been calculated'];
disp(message);
sound(sin(2*pi*10*(1:4000)/100));


