%% ULADesigner
%% Tips

%% Reset system
clear;
close all;

c = 1500;
%% Load parameters
load('START.mat','FC','RElement_Num','RSubarray_Num','TElement_Num','TSubarray_Pos',...
    'ULA_InstallationAngle','ULA_VerticalAngle');
%% ULA structure
Element_spacing = 0.5*(c/FC); % Array element spacing
ULA_Elements = RElement_Num*ones(1,RSubarray_Num);
ULA_Elements = [ULA_Elements TElement_Num];
ULA_SubarrayLength = ULA_Elements.*Element_spacing;
% Horizontal beam angle of subarray
ULA_SubarrayBeamAngle = 0.88*((1500/FC)./ULA_SubarrayLength)*180/pi;
% Array location(using the transimitting array as the Origin)
RSubarrayLength = RElement_Num*Element_spacing;
ULA_SubarrayPosition = flip((TSubarray_Pos-1)*RSubarrayLength-(RSubarray_Num-1)*RSubarrayLength:RSubarrayLength:(TSubarray_Pos-1)*RSubarrayLength);

%% Save data
save('ULADesigner.mat','ULA_SubarrayLength','ULA_SubarrayPosition','ULA_SubarrayBeamAngle',...
    'ULA_Elements','ULA_InstallationAngle','ULA_VerticalAngle');
sound(sin(2*pi*10*(1:4000)/100));

