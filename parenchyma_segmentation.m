%% SEGMENTATION
set(0,'DefaultFigureWindowStyle','normal')

%% Init, set params
clear; clc; close all;
CurrDir = pwd;

%% Data selection:
datafolder = strcat(CurrDir,'\Data');
DataFiles{1} = '\Rat_fMRI_1mm.mat';
DataFiles{2} = '\Rat_fMRI_3mm.mat';

%% Run all

load(strcat(datafolder,DataFiles{1}));

figure, imshow(coreg.MRI, []); axis on;
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
parench_mask = hFH.createMask();

for i = 1:length(DataFiles)
    save(strcat(datafolder,DataFiles{i}),'parench_mask',"-append");
end
