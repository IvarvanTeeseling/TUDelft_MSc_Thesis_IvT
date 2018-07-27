close all; clear; clc

% Read file
[Im, map] = imread('E:\CloudStation\03 Studie\MSc Thesis\09 Experiments\Run_00_Trials\02 Crack Length Processing\Images\00_02_Fatigue_Reference_30%_0_20180704_16_11_36,207.bmp');

newmap = rgb2gray(map);

% Update image
figure(1)
ImFig(1) = imshow(Im, map);

% Create Sobel-Filter
SobelF = fspecial('sobel');

% Apply Sobel filter, and extract crack by dilation, thresholding and
% opening
figure(2)
Im = imcrop(Im, map);

% Update image
figure(3)
ImFig(2) = imshow(Im);

Im = imfilter(Im, SobelF, 'replicate');

% Update image
figure(4)
ImFig(3) = imshow(Im);

elpIm = Im >100;
Im = imdilate(Im, strel('disk',2));
Im = bwareaopen(Im, 1500);

figure(5)
ImFig(4) = imshow(Im);


%[temp1,temp2] = find(Im);