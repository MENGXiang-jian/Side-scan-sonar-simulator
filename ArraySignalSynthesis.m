%% Array signal synthesis (frequency domain)
%% Tips


%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','FS','TRANSMITTER','RSUBARRAY');
load('Seabed.mat','Seabed_Z','Seabed_Zm');
load('WorkingSignal.mat','SIGNAL_fft','SIGNAL_FPnum','SIGNAL_FP','T','FR');
% Sensitivity dB(1V/uPa)
Mx = -199.4;
%% Signal synthesis process
tic
TT = T(1:round(FR/FS):end);
for transmitter = TRANSMITTER
% transmitter = 7;  
    
    filename = ['SonarParameters_',num2str(transmitter),'.mat'];
    load(filename,'PT','EL');
    filename = ['NavigationScanning_',num2str(transmitter),'.mat'];
    load(filename,'Seabed_Scanned');
    
    Signal = zeros(max(RSUBARRAY),length(TT));
    for Rsubarray = RSUBARRAY
% Rsubarray = 3;
        
        % (1) Signal amplitude
        EchoLevel = EL{1,Rsubarray};
        voltage = 10.^((EchoLevel+Mx)/20);     

        % (2) Time delay
        delay = PT{1,Rsubarray};
        
        % (3) Delay & Voltage
        delay(Seabed_Z == Seabed_Zm) = NaN; % only targets
        
        voltage(isnan(delay)) = NaN;  
        Delay = (delay(~isnan(delay)))';
        Voltage = (voltage(~isnan(voltage)))';
        
        % (4) fft
        signal_fft = zeros(1,length(SIGNAL_FP));
        for fre = 1:length(SIGNAL_FPnum)
            phase = exp(-1j*2*pi*SIGNAL_FP(fre)*Delay); 
            Phase = Voltage.*phase;
            PHASE = SIGNAL_fft(1,SIGNAL_FPnum(fre))*Phase;
            signal_fft(1,fre) = sum(PHASE(~isnan(PHASE)));
            
            message = [num2str(fre),'/',num2str(length(SIGNAL_FPnum))];
            disp(message);
        end              
        
        % (5) generation
        Signal_fft = zeros(1,length(SIGNAL_fft));
        Signal_fft(1,SIGNAL_FPnum) = signal_fft; 
        signal = real(ifft(Signal_fft)); 
        Signal(Rsubarray,:) = signal(1:round(FR/FS):end);
        
        disp(Rsubarray);
    end
    
    filename = ['ArraySignal_',num2str(transmitter),'.mat'];
    save(filename,'Signal','FR','T','FS','TT');
    disp(filename);
    sound(sin(2*pi*20*(1:4000)/100));          
end
toc

%% Save data
message = [num2str(transmitter),' groups of array signal have been synthesized'];
disp(message);
sound(sin(2*pi*10*(1:4000)/100));

