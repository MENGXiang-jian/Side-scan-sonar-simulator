%% Navigation scanning
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','NavigationSpeed_mpers','FC','DR','SP','Tracks','Sway','Heave',...
    'Roll','Yaw','Pitch','RSUBARRAY','TRANSMITTER','TSubarray_Pos','Scoring_criteria');
load('Seabed.mat','Seabed_X','Seabed_Y','Seabed_Z');
load('ULADesigner.mat','ULA_SubarrayPosition','ULA_SubarrayBeamAngle');
load('ULA_BeamPattern.mat','Angle','scale','DI_SUBARRAY');
%% Sound propagation path
for transmitter = TRANSMITTER
% % % transmitter = 31;

%     Tracks = 0*ones(1,length(TRANSMITTER)+1);
%     Sway = 0*ones(1,length(TRANSMITTER)+1);
%     Heave = 0*ones(1,length(TRANSMITTER)+1);
%     Roll = 0*ones(1,length(TRANSMITTER)+1);
%     Yaw = 0.5*ones(1,length(TRANSMITTER)+1);
%     Pitch = 25*ones(1,length(TRANSMITTER)+1);
    
    Pn_x = Sway(transmitter)+ULA_SubarrayPosition*cosd(Pitch(transmitter))*sind(Yaw(transmitter));
    Pn_y = Tracks(transmitter)+ULA_SubarrayPosition*cosd(Pitch(transmitter))*cosd(Yaw(transmitter));
    Pn_z = Heave(transmitter)+ULA_SubarrayPosition*sind(Pitch(transmitter));
    
    Pn1_x = Sway(transmitter+1)+ULA_SubarrayPosition*cosd(Pitch(transmitter+1))*sind(Yaw(transmitter+1));
    Pn1_y = Tracks(transmitter+1)+ULA_SubarrayPosition*cosd(Pitch(transmitter+1))*cosd(Yaw(transmitter+1));
    Pn1_z = Heave(transmitter+1)+ULA_SubarrayPosition*sind(Pitch(transmitter+1));    
    
    % Cache
    TD = cell(1,max(RSUBARRAY));
    RD = cell(1,max(RSUBARRAY));
    RD_Simp = cell(1,max(RSUBARRAY));
    IncidenceAngle = cell(1,max(RSUBARRAY));
    DOA = cell(1,max(RSUBARRAY));
    Seabed_Scanned = cell(1,max(RSUBARRAY));
    ENERGY = cell(1,max(RSUBARRAY));
    for Rsubarray = RSUBARRAY
% Rsubarray = 2;
        
        % (1) Propagation distance
        s1 = sqrt((Seabed_X-Pn_x(1,TSubarray_Pos)).^2+(Seabed_Y-Pn_y(1,TSubarray_Pos)).^2+(Seabed_Z-Pn_z(1,TSubarray_Pos)).^2);
        s2 = sqrt((Seabed_X-Pn_x(1,Rsubarray)).^2+(Seabed_Y-Pn_y(1,Rsubarray)).^2+(Seabed_Z-Pn_z(1,Rsubarray)).^2);
        s3 = sqrt((Seabed_X-Pn1_x(1,Rsubarray)).^2+(Seabed_Y-Pn1_y(1,Rsubarray)).^2+(Seabed_Z-Pn1_z(1,Rsubarray)).^2);
        s4 = sqrt((Pn_x(1,Rsubarray)-Pn1_x(1,Rsubarray)).^2+(Pn_y(1,Rsubarray)-Pn1_y(1,Rsubarray)).^2+(Pn_z(1,Rsubarray)-Pn1_z(1,Rsubarray)).^2);
        
        cos_alpha = (s2.^2+s4^2-s3.^2)./(2*s4.*s2);
        
        A = (c^2/NavigationSpeed_mpers^2)-1;
        B = 2*(s2.*cos_alpha-(c/NavigationSpeed_mpers)*s1);
        C = s1.^2-s2.^2;
        
        delta = (B.^2-4*A.*C);
        
        x1 = (-B+sqrt(delta))/(2*A);
        x2 = (-B-sqrt(delta))/(2*A);       
        y1 = (c/NavigationSpeed_mpers)*x1-s1;
        y2 = (c/NavigationSpeed_mpers)*x2-s1;
        
        if sum(sum(y1)) > sum(sum(y2))
            ym = y1;xm = x1;
        else
            ym = y2;xm = x2;
        end
        
        Rd = ym;         
%         Rd = s2;xm = 0; % stop&hop motion model
        
        % (2) Incidence angle
        Incidenceangle = acosd((Seabed_Z-Pn_z(Rsubarray))./Rd);
        
        % (3) Horizontal angle        
        P_x = Pn_x(Rsubarray)+(xm/s4)*(Pn1_x(Rsubarray)-Pn_x(Rsubarray));
        P_y = Pn_y(Rsubarray)+(xm/s4)*(Pn1_y(Rsubarray)-Pn_y(Rsubarray));
        P_z = Pn_z(Rsubarray)+(xm/s4)*(Pn1_z(Rsubarray)-Pn_z(Rsubarray));
        
        if Rsubarray == 5
            s4_Ref = sqrt((Pn1_x(1)-Pn_x(1))^2+(Pn1_y(1)-Pn_y(1))^2+(Pn1_z(1)-Pn_z(1))^2);
            xm_Ref = (xm/(s4/SP))*(s4_Ref/SP);
            P_x_R = Pn_x(1)+(xm_Ref/s4_Ref)*(Pn1_x(1)-Pn_x(1));
            P_y_R = Pn_y(1)+(xm_Ref/s4_Ref)*(Pn1_y(1)-Pn_y(1));
            P_z_R = Pn_z(1)+(xm_Ref/s4_Ref)*(Pn1_z(1)-Pn_z(1));
        else
            s4_Ref = sqrt((Pn1_x(5)-Pn_x(5))^2+(Pn1_y(5)-Pn_y(5))^2+(Pn1_z(5)-Pn_z(5))^2);
            xm_Ref = (xm/(s4/SP))*(s4_Ref/SP);
            P_x_R = Pn_x(5)+(xm_Ref/s4_Ref)*(Pn1_x(5)-Pn_x(5));
            P_y_R = Pn_y(5)+(xm_Ref/s4_Ref)*(Pn1_y(5)-Pn_y(5));
            P_z_R = Pn_z(5)+(xm_Ref/s4_Ref)*(Pn1_z(5)-Pn_z(5));            
        end
        
        ym_Ref = sqrt((Seabed_X-P_x_R).^2+(Seabed_Y-P_y_R).^2+(Seabed_Z-P_z_R).^2);
        s5 = sqrt((P_x_R-P_x).^2+(P_y_R-P_y).^2+(P_z_R-P_z).^2);
        
        Doa = abs(90-(acosd((s5.^2+Rd.^2-ym_Ref.^2)./(2.*Rd.*s5))));
        
        % (4) Modeling simplified algorithm (Modified energy function based on acoustic principle )       
        % ① Shadow weight coefficients
        w_Shadow = ones(size(Seabed_Z));
        Theat = 90-(acosd((s5.^2+Rd.^2-ym_Ref.^2)./(2.*Rd.*s5)));
        delt_doa = 0.05;
        ID = zeros(1,2);
        NUM = reshape(1:numel(Theat),size(Theat,1),size(Theat,2));
        for doa = min(min(Theat)):delt_doa:max(max(Theat)) % min(min(Theat))
            THEAT = Theat;THEAT(THEAT<doa) = NaN;THEAT(THEAT>doa+delt_doa) = NaN;  
            
            Location = sqrt((Seabed_X-Pn_x(1,TSubarray_Pos)).^2+(Seabed_Y-Pn_y(1,TSubarray_Pos)).^2); Location(isnan(THEAT)) = NaN;
            location = Location(~isnan(Location));
            [~,id] = ismember(unique(location),Location);
            
            num = NUM(id);           
            incidenceangle = Incidenceangle(id);
            if length(incidenceangle) > 3
                [value,Index] = findpeaks(incidenceangle);            
                M = zeros(1,2);
                for ii = 1:length(Index)
                    index = find(incidenceangle(Index(ii):end) < value(ii)); M = [M,(Index(ii)+index-1)'];
                end
                MM = unique(M);MM = MM(MM ~= 0);           
                [~,Id] = ismember(num(MM),NUM); ID = [ID,Id'];
            end
            disp(doa);
        end
        
        ID = unique(ID);
        w_Shadow(ID(ID ~= 0)) = NaN;        
%         Seabed_Z(ID(ID ~= 0)) = NaN;
        
        % ② Scattering energy component
        E_Scattering = abs(reshape(mapminmax(Incidenceangle(:)',-1,0),size(Incidenceangle)));
        % ③ Curvature energy component
        [Gaussian,Average,Pmax,Pmin] = surfature(Seabed_X,Seabed_Y,Seabed_Z);
        E_Curvature = mapminmax(abs(Pmax),0,1);
%         E_Curvature = reshape(mapminmax(Average_abs(:)',0,1),size(Average_abs));
        % ④ View energy component     
        E_View = ones(size(Seabed_Z));
        View = DI_SUBARRAY(round(Doa/scale)+9001);           
        E_View = reshape(mapminmax(View(:)',0,1),size(View)); % 
        E_View(Rd > DR) = NaN;
        E_View(Doa > 0.5*ULA_SubarrayBeamAngle(1,Rsubarray)) = NaN;      
        % ⑤ Modified energy function
        Energy = w_Shadow.*(1*E_Scattering+1*E_Curvature+4*E_View);
        
        Energy_Score = reshape(mapminmax(Energy(:)',0,1),size(Energy));
        Seabed_scanned = ones(size(Seabed_Z));
        Seabed_scanned(isnan(Energy_Score)) = NaN;
        Seabed_scanned(Energy_Score < Scoring_criteria) = NaN;        
%         Energy_Score(Energy_Score < Scoring_criteria) = NaN;       
        
        % (5) Data
        TD{1,Rsubarray} = s1.*Seabed_scanned;
        RD{1,Rsubarray} = Rd.*Seabed_scanned;
        IncidenceAngle{1,Rsubarray} = Incidenceangle.*Seabed_scanned;
        DOA{1,Rsubarray} = Doa.*Seabed_scanned;
        Seabed_Scanned{1,Rsubarray} = Seabed_scanned; 
        ENERGY{1,Rsubarray} = Energy_Score; 
        
        % (6) Display
% % %         figure(1);        
% % %         Display = Seabed_Z.*Seabed_scanned;
% % %         mesh(Seabed_X(1,:),Seabed_Y(:,1),Display);
% % %         shading interp;
% % %         xlabel('Transverse direction（m）','FontSize',15); 
% % %         ylabel('Navigation direction（m）','FontSize',15);
% % %         zlabel('Incidence angle（°）','FontSize',15);
% % %         set(gca,'ZDir','reverse','FontSize',15);
% % %         zlim([0,30]);  
% % %         
% % %         disp(Rsubarray);
    end
    
    % Save data
    filename = ['NavigationScanning_',num2str(transmitter),'.mat'];
    save(filename,'TD','RD','IncidenceAngle','DOA','Seabed_Scanned','ENERGY'); 
    disp(filename);
    sound(sin(2*pi*20*(1:4000)/100));
end

%% Save data
message = [num2str(transmitter),' groups of sound propagation parameters have been calculated'];
disp(message);
sound(sin(2*pi*10*(1:4000)/100));


