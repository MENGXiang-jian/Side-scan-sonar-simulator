%% Analog circuit processing
%% Tips


%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','noise_fl','noise_fh','Fixed_gain_dB','FC','TRANSMITTER','RSUBARRAY');
%% Gain processing
% ① TVG
load('ArraySignal_1.mat','FS','TT');
alpha=((0.1*(FC/1000)^2)/(1+(FC/1000)^2)+(40*(FC/1000)^2)/(4100+(FC/1000)^2)...
    +2.75*(10^-4)*(FC/1000)^2+0.003)*1;
tvg_distance = TT*c/2;
tvg = 2*(20*log10(1.093613298*tvg_distance)+alpha*((1.093613298*tvg_distance)/1000));
TVG_gain = 10.^(tvg/20);
TVG_gain(1) = 0;
% ② Fixed Gain
Fixed_gain = 10^(Fixed_gain_dB/20);
%% Noise
noise = real(exp(1j*2*pi*noise_fl*TT))+real(exp(1j*2*pi*noise_fh*TT));

%% Receiver simulation

for transmitter = TRANSMITTER
%     transmitter = 31;
    
    filename = ['ArraySignal_',num2str(transmitter),'.mat'];
    load(filename,'Signal','FR','T','FS','TT'); 
    filename = ['SonarParameters_',num2str(transmitter),'.mat'];
    load(filename,'DT');
    
    % Cache
    Asignal_Original = zeros(length(RSUBARRAY),length(TT));
    Asignal_Transducer = zeros(length(RSUBARRAY),length(TT));
    Asignal_BandPassFilter = zeros(length(RSUBARRAY),length(TT));
    Asignal_TVG = zeros(length(RSUBARRAY),length(TT));
    Asignal = zeros(length(RSUBARRAY),length(TT));
    for Rsubarray = RSUBARRAY
% Rsubarray = 2;       
        % 0) Original signal
        Asignal_Original(Rsubarray,:) = Signal(Rsubarray,:);
        % 1) Receiving signal
        % Noise        
        Amplitude_max = max(Asignal_Original(Rsubarray,:));                
        Noise = 0*(Amplitude_max/(10^(DT/10)))*(noise)...
        +0*awgn(Asignal_Original(Rsubarray,:),66,'measured')...
        +1*(Amplitude_max/(10^(DT/20)))*randn(1,length(Asignal_Original(Rsubarray,:)));   
        Noise = filter(AChebyshevTypeII_BandPass,Noise);
        Asignal_Transducer(Rsubarray,:) = 1*Asignal_Original(Rsubarray,:)+0*Noise;
        % 2) Bandpass filtering
        Asignal_BandPassFilter(Rsubarray,:) = filter(AChebyshevTypeII_BandPass,Asignal_Transducer(Rsubarray,:));
        % 3) TVG processing
        Asignal_TVG(Rsubarray,:)  = Asignal_BandPassFilter(Rsubarray,:).*TVG_gain;     
        % 4) Fixed gain processing
        Asignal(Rsubarray,:) = Asignal_TVG(Rsubarray,:)*Fixed_gain;        
        
        disp(Rsubarray);
    end
    
    filename = ['AnalogCircuitProcessing_',num2str(transmitter),'.mat']; % transmitter
    save(filename,'Asignal_Original','Asignal_Transducer','Asignal_BandPassFilter','Asignal_TVG','Asignal','T','FR','TT','FS','Noise','TVG_gain','Fixed_gain');
    disp(filename);
    sound(sin(2*pi*20*(1:4000)/100));       
end

%% Save data
message = [num2str(transmitter),' groups of analog processing have been synthesized'];
disp(message);
sound(sin(2*pi*10*(1:4000)/100));

