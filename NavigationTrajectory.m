%% Motion error setting
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','NavigationSpeed_mpers','Tracks','Sway','Heave','Roll','Yaw','Pitch');

%% Display
figure(1)
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);

subplot(2,1,1);
plot(Tracks,Tracks);
hold on;
plot(Tracks,Sway);
hold on;
plot(Tracks,Heave);
set(gca,'FontSize',14);
legend('Tracks','Sway','Heave')
title('Displacement');
hold off;

subplot(2,1,2);
plot(Tracks,Yaw);
hold on;
plot(Tracks,Pitch);
hold on;
plot(Tracks,Roll);
set(gca,'FontSize',14);
legend('Yaw','Pitch','Roll')
title('Attitude');
hold off;

figure(2)
scrsz = [20,40,1500,700];
set(gcf,'Position',scrsz);

load('START.mat','Tracks','TSubarray_Pos');
load('ULADesigner.mat','ULA_SubarrayLength','ULA_SubarrayPosition','ULA_SubarrayBeamAngle');
Pn_x = zeros(1,length(Tracks)-1);
Pn_y = zeros(1,length(Tracks)-1);
Pn_z = zeros(1,length(Tracks)-1);
for transmitter = 1:length(Tracks)-1
    Pn_x(transmitter) = Sway(transmitter)+(ULA_SubarrayPosition(TSubarray_Pos).*sind(Yaw(transmitter)));
    Pn_y(transmitter) = Tracks(transmitter)+(ULA_SubarrayPosition(TSubarray_Pos).*cosd(Yaw(transmitter)).*cosd(Pitch(transmitter)));
    Pn_z(transmitter) = Heave(transmitter)+ULA_SubarrayPosition(TSubarray_Pos).*sind(Pitch(transmitter));
end

plot3(Pn_x,Pn_y,Pn_z,'o-');
grid on;
view(60,45);
xlim([-5,5]);  
% axis equal;
%% Save data

sound(sin(2*pi*10*(1:4000)/100));


