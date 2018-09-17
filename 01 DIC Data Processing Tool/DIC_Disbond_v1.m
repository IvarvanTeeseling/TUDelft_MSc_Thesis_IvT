%% Initialize MATLAB
clear; close all; clc

% Set timer
tic

%% File selection - 2_11_A_R_26

% Path to DIC image files
dicPathSource = fullfile('C:\Users\Ivar\Desktop\tmp');
%dicPathSource = fullfile('D:\Thesis Experimental Data\Run_02_CLS\2_12_A_HL_26\DIC\Matlab Data');


% First loaded test image
dicFirstLoadImage = 5;
camFirstLoadImage = 1;

% First test image for plotting
testFirstPlotImage = 15;

% # Load cycles at the first, loaded image
nFirstImage = 500;

% # Load cycles / image
nStep = 500;

% Rectangle coordinates
b = 44;
lA = 145;
StepSizeGap = 1;
rec = [-b/2 0 b lA];

% Convalution filter settings
Ksize1D = 5;

% Hi-Lock dimensions and position
Rc = 9/2+1;
xc = 0;
yc = 36;

% Frame # w.r.t. dicFirstLoadImage to be subtracted from all DIC results
% Available settings:
%   SubIm = NaN : switched off
SubIm = NaN;

%% Processing - DO NOT CHANGE

% Get all files
dicMatData = dir(fullfile(dicPathSource, '*.mat'));
dicFrameIndex = 1:length(dicMatData);

% Initiate a waitbar
f = waitbar(1/length(dicFrameIndex), 'Uploading MATLAB data. Please wait...');

% Read, extract and store the DIC .mat files selected for processing
for i = 1:length(dicFrameIndex)
    
    % DIC image index
    dicInd = dicFrameIndex(i);
    
    % Load DIC data
    dic_rw(i) = load(fullfile(dicMatData(dicInd).folder, dicMatData(dicInd).name));
    
    % Pre-allocate memory
    if i == 1
        X = zeros([size(dic_rw(i).X) length(dicFrameIndex)]);
        Y = zeros([size(dic_rw(i).X) length(dicFrameIndex)]);
        ey = zeros([size(dic_rw(i).X) length(dicFrameIndex)]);
    end
    
    % Store data in matrix form for easier handeling
    X(:,:,i) = dic_rw(i).X;
    Y(:,:,i) = dic_rw(i).Y;
    ey(:,:,i) = dic_rw(i).eyy;
    
    % Update waitbar
    waitbar(i/length(dicFrameIndex), f);
end

% Close waitbar
close(f);

%% Plot data - All DIC eyy contour plots & Video

% Select the variable to be plotted
var = ey;

% Loop through every image and store into the video
for i = 50
    
    % Extract data to be plotted
    if isnan(SubIm)
        tmpZ = var(:,:,i);
    else
        tmpZ = var(:,:,i+1)-var(:,:,1);
    end
    tmpX = X(:,:,i+(~isnan(SubIm)));
    tmpY = Y(:,:,i+(~isnan(SubIm)));
    
    % Adjust (X, Y) coordinates
    tmpX = tmpX-(min(tmpX(:))-rec(1))+StepSizeGap;
    tmpY = tmpY-(min(tmpY(:))-rec(2))+StepSizeGap;
    
    % Create interpolated xq, yq and zy values
    xc = Rc*cos(linspace(0, 2*pi, 200))+xc;
    yc = Rc*sin(linspace(0, 2*pi, 200))+yc;
    [xq, yq, zq] = dicgriddata(tmpX, tmpY, tmpZ, 1, 1, 'nofill', xc, yc);
    
    % Apply concolution filter if requested
    zq(zq==0) = NaN;
    zqsmooth = nanconv(zq, fspecial('average', 9), 'edge');

    % Find the z derrivatives
    xq_di = xq(2:end,:);
    yq_di = yq(2:end,:);
    xq_dii = xq(3:end,:);
    yq_dii = yq(3:end,:);
    zq_di = diff(zq, 1);
    zq_dii = diff(zq_di, 1); 
    zq_di = smoothdata(zq_di, 1, 'movmean', Ksize1D, 'omitnan');
    zq_dii = smoothdata(zq_dii, 1, 'movmean', Ksize1D, 'includenan');
    
    % Plot 1D curve
    q = 1;
    
    % z peaks (minima and maxima)
    [maxi, locsmaxi] = findpeaks(zq(:,q));
    [mini, locsmini] = findpeaks(zq(:,q).^-1);
    mm = [maxi locsmaxi ones(size(maxi)) ; ...
        mini.^-1 locsmini zeros(size(mini))];
    mm_sort = sortrows(mm, 2);
    
    mm_diff = diff(mm_sort(1+mm_sort(1,3):end,1));
    
    locind = mm_sort(find(mm_diff == max(mm_diff))+mm_sort(1,3), 2);
    
    locyq = yq(locind, q);
    loczq = zq(locind, q);
    
    h =  findobj('type','figure');
    n = length(h);
    figure(n+1)
    hold on
    plot(yq(:,q), zq(:,q))
    plot(yq(:,q), zqsmooth(:,q))
    plot(yq(mm(:,2), q), mm(:,1), 'o')
    plot(locyq, loczq, 'rs')
    legend('z_{interp.}', 'z_{interp. + mean filter}', 'disbond location','z_{peaks}')
    hold off
    
    [zb, bA] = locateDisbondFront(yq, zqsmooth, 'findvalley');
    
    
    figure(2)
    subplot(1,3,1)
    contourf(xq, yq, zq)
    subplot(1,3,2)
    contourf(xq, yq, zqsmooth)
    subplot(1,3,3)
    contourf(xq, yq, zb)
    
    % Perform z thresholding
    
    zq_th = zq;
    zq_th(zq_th==min(zq_th, [], 1)) = 1;
    zq_th(zq_th~=1) = 0;
    
    % zdii thresholding
    zdii_th = zq_dii;
    zdii_th(zdii_th==max(zdii_th, [], 1)) = 1;
    zdii_th(zdii_th~=1) = 0;
    
    % Plot the contour data
    tmpZ(tmpZ==0) = NaN;
    tmpX(isnan(tmpZ)) = NaN;
    tmpY(isnan(tmpZ)) = NaN;
    
    h =  findobj('type','figure');
    n = length(h);
    figure(n+1)
    subplot(1,7,1)
    contourf(tmpX, tmpY, tmpZ)
    title('Raw')
    subplot(1,7,2)
    contourf(xq, yq, zq)
    title('Interp.')
    subplot(1,7,3)
    contourf(xq, yq, zqsmooth)
    title('Interp. + filter')
    subplot(1,7,4)
    contourf(xq_di, yq_di, zq_di)
    title('(Interp. + filter) dy/dx')
    subplot(1,7,5)
    contourf(xq_dii, yq_dii, zq_dii)
    title('(Interp. + filter) d^2y/dx^2')
    subplot(1,7,6)
    contourf(xq_di, yq_di, zdi_th)
    title('Disbond dy/dx')
    subplot(1,7,7)
    contourf(xq_dii, yq_dii, zdii_th)
    title('Disbond d^2y/dx^2')
end
