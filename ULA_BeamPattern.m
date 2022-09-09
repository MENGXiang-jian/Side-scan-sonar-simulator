%% ULA_BeamPattern
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
% % % load('ULA_DesignParameters.mat','ULA_Elements','ULA_Length','ULA_BeamAngle');
% % % load('Signal.mat','Chirp_frequency','ULA_SlantDistance');
load('START.mat','FC','FR','BW','PW','DR','SP');
load('ULADesigner.mat','ULA_SubarrayLength','ULA_SubarrayPosition','ULA_SubarrayBeamAngle',...
    'ULA_Elements','ULA_InstallationAngle','ULA_VerticalAngle');
ULA_WholeRLength = sum(ULA_SubarrayLength(1,1:end-1));
lamda = c/FC;

%% Spherical wave model
% Fresnel-Fraunhofer
Range_Fresnel = (ULA_WholeRLength^2)/lamda;
Range_Fraunhofer = (pi*ULA_WholeRLength^2)/lamda;

% Angle
scale = 0.01;
Angle = -90:scale:90; 
% Distance
Range = 1.356*ULA_WholeRLength:5:DR;
Range = 2500;

%% Element beam pattern
L_element = lamda/2;

DI_element = zeros(1,length(Angle));
for aa = 1:length(Angle)     
    DI_element(1,aa) = abs(sin(pi*FC*L_element/c*sind(Angle(aa)))/(pi*FC*L_element/c*sind(Angle(aa))));
end
DI_element((length(Angle)+1)/2) = 1;

%% ULA beam pattern

DI_elementsSumR = zeros(length(Range),length(Angle));
DI_SubarrayR = zeros(length(Range),length(Angle));
DI_SUBARRAY = zeros(length(Range),length(Angle));

DI_arraySum = zeros(length(Range),length(Angle));
DI_Array = zeros(length(Range),length(Angle));
DI_ARRAYr = zeros(length(Range),length(Angle));

DI_elementsSumT = zeros(length(Range),length(Angle));
DI_SubarrayT = zeros(length(Range),length(Angle));
DI_ARRAYt = zeros(length(Range),length(Angle));

for rr = 1:length(Range)

    for aa = 1:length(Angle)

       %% Receiving subarray
        Dr = abs((floor(ULA_Elements(1)/2)+1)-(1:ULA_Elements(1)))*L_element;
        PropagationDelay_elementsR(1:floor(ULA_Elements(1)/2)) = (sqrt(Range(rr)^2+Dr(1:floor(ULA_Elements(1)/2)).^2-2*Range(rr).*Dr(1:floor(ULA_Elements(1)/2))*cosd(90+Angle(aa))))/c;
        PropagationDelay_elementsR(floor(ULA_Elements(1)/2)+1:ULA_Elements(1)) = (sqrt(Range(rr)^2+Dr(floor(ULA_Elements(1)/2)+1:ULA_Elements(1)).^2-2*Range(rr).*Dr(floor(ULA_Elements(1)/2)+1:ULA_Elements(1))*cosd(90-Angle(aa))))/c;
        
        DI_elementsR = exp(1j*(2*pi*FC*PropagationDelay_elementsR));
        DI_elementsSumR(rr,aa) = sum(real(DI_elementsR));
        DI_SubarrayR(rr,aa) = DI_elementsSumR(rr,aa).*DI_element(aa);
        
       %% Receiving array
        DR = abs((floor((length(ULA_SubarrayLength)-1)/2)+1)-(1:(length(ULA_SubarrayLength)-1)))*ULA_SubarrayLength(1);
        PropagationDelay_Subarray(1:floor((length(ULA_SubarrayLength)-1)/2)) = (sqrt(Range(rr)^2+DR((1:floor((length(ULA_SubarrayLength)-1)/2))).^2-2*Range(rr).*DR((1:floor((length(ULA_SubarrayLength)-1)/2)))*cosd(90+Angle(aa))))/c;
        PropagationDelay_Subarray(floor((length(ULA_SubarrayLength)-1)/2)+1:(length(ULA_SubarrayLength)-1)) = (sqrt(Range(rr)^2+DR(floor((length(ULA_SubarrayLength)-1)/2)+1:(length(ULA_SubarrayLength)-1)).^2-2*Range(rr).*DR(floor((length(ULA_SubarrayLength)-1)/2)+1:(length(ULA_SubarrayLength)-1))*cosd(90-Angle(aa))))/c;
        
        DI_array = exp(1j*(2*pi*FC*PropagationDelay_Subarray)); 
        DI_arraySum(rr,aa) = sum(real(DI_array));     
        DI_Array(rr,aa) = DI_arraySum(rr,aa).*DI_SubarrayR(rr,aa);             
        
        %% Transmitting array
        DT = abs((floor(ULA_Elements(end)/2)+1)-(1:ULA_Elements(end)))*L_element;
        PropagationDelay_elementsT(1:floor(ULA_Elements(end)/2)) = (sqrt(Range(rr)^2+DT(1:floor(ULA_Elements(end)/2)).^2-2*Range(rr).*DT(1:floor(ULA_Elements(end)/2))*cosd(90+Angle(aa))))/c;
        PropagationDelay_elementsT(floor(ULA_Elements(end)/2)+1:ULA_Elements(end)) = (sqrt(Range(rr)^2+DT(floor(ULA_Elements(end)/2)+1:ULA_Elements(end)).^2-2*Range(rr).*DT(floor(ULA_Elements(end)/2)+1:ULA_Elements(end))*cosd(90-Angle(aa))))/c;
        
        DI_elementsT = exp(1j*(2*pi*FC*PropagationDelay_elementsT));
        DI_elementsSumT(rr,aa) = sum(real(DI_elementsT));        
        DI_SubarrayT(rr,aa) = DI_elementsSumT(rr,aa).*DI_element(aa);          

    end
    
    %% Normalization processing
    DI_ELE = 20*log10(abs(DI_element)/max(abs(DI_element)));
    DI_SUBARRAY(rr,:) = 20*log10(abs(DI_SubarrayR(rr,:))/max(abs(DI_SubarrayR(rr,:))));   
    DI_ARRAYr(rr,:) = 20*log10(abs(DI_Array(rr,:))/max(abs(DI_Array(rr,:))));                                        
    DI_ARRAYt(rr,:) = 20*log10(abs(DI_SubarrayT(rr,:))/max(abs(DI_SubarrayT(rr,:))));
    
end

%% Save data
save('ULA_BeamPattern.mat','Range','Angle','scale','DI_ELE','DI_SUBARRAY','DI_ARRAYr','DI_ARRAYt');

%% Display
rr = 1;

%% (1) Element
figure(1);
% scrsz = [20,40,1500,700];
% set(gcf,'Position',scrsz); 
plot(Angle,DI_ELE,'b-'); 
xlabel('Angle(°)','fontsize',14);
ylabel('Directional sensitivity(dB)','fontsize',14);
set(gca,'FontSize',14);
% title('单个条状阵元的指向性');
% ylim([-40 0]);
% xlim([-0.25 0.25]);

%% (2) Receiving subarray
figure(2);
% scrsz = [20,40,1500,700];
% set(gcf,'Position',scrsz); 
plot(Angle,DI_SUBARRAY(rr,:),'b-'); 
xlabel('Angle(°)','fontsize',14);
ylabel('Directional sensitivity(dB)','fontsize',14);
set(gca,'FontSize',14);
% title('单个接收子阵的指向性');
% ylim([-40 0]); % 水平开角 0.5°
% xlim([-0.25 0.25]);

%% (3) Receiving array
figure(3);
% scrsz = [20,40,1500,700];
% set(gcf,'Position',scrsz); 
plot(Angle,DI_ARRAYr(rr,:),'b-'); 
xlabel('Angle(°)','fontsize',14);
ylabel('Directional sensitivity(dB)','fontsize',14);
set(gca,'FontSize',14);
% title('完整接收阵列的指向性');
% ylim([-40 0]); % 水平开角 0.5°
% xlim([-0.25 0.25]);

%% (4) Transmitting array
figure(4);
% scrsz = [20,40,1500,700];
% set(gcf,'Position',scrsz); 
plot(Angle,DI_ARRAYt(rr,:),'b-'); 
xlabel('Angle(°)','fontsize',14);
ylabel('Directional sensitivity(dB)','fontsize',14);
set(gca,'FontSize',14);
% title('单个发射阵列的指向性');
% ylim([-40 0]); % 水平开角 0.5°
% xlim([-0.25 0.25]);




