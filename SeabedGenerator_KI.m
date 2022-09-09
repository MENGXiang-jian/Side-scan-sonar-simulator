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

%% Seabed target
CURVE = -real(25*exp(1j*2*pi*0.3*Seabed_Y(:,1)));
curve_X = [100 200];
for yy = 1:size(Seabed_Y,1)    
    for curve = 1:length(curve_X)
        xx = round(curve_X(curve)/Seabed_delt_xm)+1+round(CURVE(yy)/Seabed_delt_xm);
        zz = Seabed_Zm-5;
        Seabed_Z(yy,xx:xx+10/Seabed_delt_xm) = zz;  
    end
end

% Line
Line_X = [140 240];
Seabed_Z(:,Line_X(1)/Seabed_delt_xm+1:Line_X(1)/Seabed_delt_xm+1+5/Seabed_delt_xm) = Seabed_Zm-1;
Seabed_Z(:,Line_X(2)/Seabed_delt_xm+1:Line_X(2)/Seabed_delt_xm+1+5/Seabed_delt_xm) = Seabed_Zm-1;

% % % % % Random point
% % % % num = 0;aa = numel(Seabed_Z);
% % % % Seabed_Z(datasample((find(Seabed_Z == Seabed_Zm))',num)) = Seabed_Zm-rand(1,num)*2;

%% TM generator
TM_Original = Seabed_Z(:,70/Seabed_delt_xm+1:150/Seabed_delt_xm+1);
[TM_X,TM_Y] = meshgrid(0:Seabed_delt_xm:80,0:Seabed_delt_ym:Seabed_Ym);
% surf2stl('TM_Original.stl',TM_X,TM_Y,TM_Original); % generate stl file

%% Save data
save('Seabed.mat','Seabed_Zm','Seabed_X','Seabed_Y','Seabed_Z','Seabed_delt_xm','Seabed_delt_ym',...
    'TM_X','TM_Y','TM_Original');
sound(sin(2*pi*10*(1:4000)/100));

%% Display
figure(1);

surf(Seabed_X(1,:),Seabed_Y(:,1),Seabed_Z);
shading interp;
xlabel('Transverse direction（m）','FontSize',16); 
ylabel('Navigation direction（m）','FontSize',16);
zlabel('Water depth（m）','FontSize',16);
set(gca,'ZDir','reverse','FontSize',16);
zlim([0,40]);  
% axis equal;

figure(2);

surf(TM_X(1,:),TM_Y(:,1),TM_Original);
shading interp;
xlabel('Transverse direction（m）','FontSize',16); 
ylabel('Navigation direction（m）','FontSize',16);
zlabel('Water depth（m）','FontSize',16);
set(gca,'ZDir','reverse','FontSize',16);
zlim([0,40]);  
% axis equal;






