%% Seabed model
%% Tips

%% Reset system
clear;
close all;

%% Load point cloud
load('Seabed.mat','Seabed_Zm','Seabed_X','Seabed_Y','Seabed_delt_xm','Seabed_delt_ym',...
    'TM_X','TM_Y','TM_Original');

load('TM_xyz.txt'); % point cloud 
x = round(TM_xyz(:,1)/Seabed_delt_xm)+1;
y = round(TM_xyz(:,2)/Seabed_delt_ym)+1;
z = TM_xyz(:,3);

Seabed_Z = Seabed_Zm*ones*(size(Seabed_X));
for num = 1:length(z)    
    Seabed_Z(y(num),x(num)) = z(num);    
end

save('Seabed.mat','Seabed_Zm','Seabed_X','Seabed_Y','Seabed_Z','Seabed_delt_xm','Seabed_delt_ym',...
    'TM_X','TM_Y','TM_Original');
%% Load TM
% [tri,fileform,A,S] = stlread('TM.stl');

fig = figure(1);
model = stl2matlab('TM.stl');
patch(model{1},model{2},model{3},'b');
xlabel('Transverse direction（m）','FontSize',16); 
ylabel('Navigation direction（m）','FontSize',16);
zlabel('Water depth（m）','FontSize',16);
set(gca,'ZDir','reverse','FontSize',16);
% zlim([0,40]); 
% axis equal; 
%% Display
figure(2);
surf(TM_X(1,:),TM_Y(:,1),Seabed_Z);
shading interp;
xlabel('Transverse direction（m）','FontSize',16); 
ylabel('Navigation direction（m）','FontSize',16);
zlabel('Water depth（m）','FontSize',16);
set(gca,'ZDir','reverse','FontSize',16);
zlim([0,40]);  
% axis equal;




