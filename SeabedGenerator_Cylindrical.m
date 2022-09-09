%% Seabed model
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameter
load('START.mat','NavigationSpeed_mpers','Voyage');
%% Build model
Seabed_Xm = 260;
Seabed_Ym = Voyage;
Seabed_Zm = 30;

Seabed_delt_xm = 0.1; 
Seabed_delt_ym = 0.01; 

Seabed_x = round(Seabed_Xm/Seabed_delt_xm)+1;
Seabed_y = round(Seabed_Ym/Seabed_delt_ym)+1;

Seabed_Z = Seabed_Zm*ones(Seabed_y,Seabed_x);
[Seabed_X,Seabed_Y] = meshgrid(0:Seabed_delt_xm:Seabed_Xm,0:Seabed_delt_ym:Seabed_Ym);

%% Semi-cylindrical target
SC_Ym = [2 5 8];
SC_Xm = [30 90 150 210];

for yy = 1:length(SC_Ym)
    for xx = 1:length(SC_Xm)        
        % target position
        SC_ym = SC_Ym(yy);
        SC_xm = SC_Xm(xx);
        % target size
        R_SC = 5;
        L_SC = 1;
        % build target
        SC_x = 0:Seabed_delt_xm:2*R_SC;
        SC_y = 0:Seabed_delt_ym:L_SC;
        [SC_x,SC_y] = meshgrid(SC_x,SC_y); 
        SC_Zm = real(sqrt(R_SC^2-(abs(R_SC-SC_x)).^2)); 
        SC_Zm = Seabed_Zm-SC_Zm;        
        % Place target
        SC_xs = round((SC_xm-R_SC)/Seabed_delt_xm)+1;
        SC_xe = round((SC_xm+R_SC)/Seabed_delt_xm)+1; 
        SC_ys = round((SC_ym-0.5*L_SC)/Seabed_delt_ym)+1; 
        SC_ye = round((SC_ym+0.5*L_SC)/Seabed_delt_ym)+1;   
        
        Seabed_Z(SC_ys:SC_ye,SC_xs:SC_xe) = SC_Zm;     
    end
end

%% Save data
save('Seabed.mat','Seabed_Zm','Seabed_X','Seabed_Y','Seabed_Z','Seabed_delt_xm','Seabed_delt_ym');
sound(sin(2*pi*10*(1:4000)/100));

%% Display
figure(1);
% scrsz = [20,40,1500,700];
% set(gcf,'Position',scrsz);

surf(Seabed_X(1,:),Seabed_Y(:,1),Seabed_Z);
shading interp;
xlabel('Transverse direction（m）','FontSize',16); 
ylabel('Navigation direction（m）','FontSize',16);
zlabel('Water depth（m）','FontSize',16);
set(gca,'ZDir','reverse','FontSize',16);
zlim([0,40]);  
% axis equal;








