%% Initialize MATLAB
clear; close all; clc

% Set timer
tic

%% File selection

% CLS specimen ID
%   1 = 01_03_Fat_Ref_30kN > No taper
%   2 = 2_11_A_R_26
%   3 = 2_12_A_HL_26
%   4 = 2_13_A_HL_26
clsID = 3;

% Path to camera image files
camPathSource = {
    'D:\Thesis Experimental Data\Run_01_CLS\01_03_Fat_Ref_30kN\Data Processing' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_11_A_R_26\Data Processing' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_12_A_HL_26\Data Processing' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_13_A_HL_26\Data Processing'
    };

camNameL = {
    'Camera Post Processing - 01_03_Fat_Ref_30kN (Right Camera).txt' ; ...
    'Camera Post Processing - 2_11_A_R_26 (Left Camera).txt' ; ...
    'Camera Post Processing - 2_12_A_HL_26 (Left Camera).txt' ; ...
    'Camera Post Processing - 2_13_A_HL_26 (Left Camera).txt'
    };

camNameR = {
    'Camera Post Processing - 01_03_Fat_Ref_30kN (Left Camera).txt' ; ...
    'Camera Post Processing - 2_11_A_R_26 (Rigth Camera).txt' ; ...
    'Camera Post Processing - 2_12_A_HL_26 (Right Camera).txt'
    'Camera Post Processing - 2_13_A_HL_26 (Right Camera).txt'
    };

dicPathSource = {
    'D:\Thesis Experimental Data\Run_01_CLS\01_03_Fat_Ref_30kN\DIC\MATLAB Data' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_11_A_R_26\DIC\Matlab Data' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_12_A_HL_26\DIC\MATLAB Data' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_13_A_HL_26\DIC\MATLAB Data\Complete'
    };

pathTarget = {
    'D:\Thesis Experimental Data\Run_01_CLS\01_03_Fat_Ref_30kN\Data Processing' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_11_A_R_26\Data Processing' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_12_A_HL_26\Data Processing' ; ...
    'D:\Thesis Experimental Data\Run_02_CLS\2_13_A_HL_26\Data Processing'
    };

clsName = {
    '01_03_Fat_Ref_30kN' ; ...
    '2_11_A_R_26' ; ...
    '2_12_A_HL_26' ; ...
    '2_13_A_HL_26'
    };

% First loaded test image
dicLoaded = [1 5 6 5];
camLoaded = [1 1 1 1];

% First test image for plotting
testFirst = [1 1 1 5] ;

% # Load cycles at the first, loaded image
FirstImage = [500 500 500];

% # Load cycles / image
nStep = [500 500 500 500];

% Rectangle coordinates
dicEdgeGap = [1 1 1 1];
rec = [-48/2 0 48 145 ; ...
    -44/2 0 44 145 ; ...
    -44/2 0 44 145 ; ...
    -44/2 0 44 145 ; ...
    ];

circ = [0 0 0 ; ...
    0 0 0 ; ...
    10/2+1 0 35 ; ...
    10/2+1 0 45
    ];

n = [1 1 ; ...
    1 1 ; ...
    2 2 ; ...
    1 1
    ];

Ksize1D = [3 ; ...
    10 ; ...
    10 ; ...
    10
    ];

% Frame # w.r.t. dicLoaded to be subtracted from all DIC results
% Available settings:
%   SubIm = NaN : switched off
SubIm = [NaN NaN NaN -1];


%% Initialize settings

camPathSource       = camPathSource{clsID}
camNameL            = camNameL{clsID}
camNameR            = camNameR{clsID}
dicPathSource       = dicPathSource{clsID}
pathTarget          = pathTarget{clsID}
clsName             = clsName{clsID}
dicLoaded           = dicLoaded(clsID)
camLoaded           = camLoaded(clsID)
testFirst           = testFirst(clsID)
FirstImage          = FirstImage(clsID)
nStep               = nStep(clsID)
dicEdgeGap          = dicEdgeGap(clsID)
rec                 = rec(clsID, :)
SubIm               = SubIm(clsID)
circ                = circ(clsID,:)
n                   = n(clsID,:)
Ksize1D             = Ksize1D(clsID)

%% Upload results

% Extract DIC matlab data
dicFiles = dir(fullfile(dicPathSource, '*.mat'));

% Extract camera image disbond length data
camDataL = textread(fullfile(camPathSource, camNameL));
camDataR = textread(fullfile(camPathSource, camNameR));

% Only keep measurements (images) by both left and right camera
camFrIndex = intersect(camDataL(:,2)', camDataR(:,2)')';

% Ingore camera images not requested for processing
camFirst = camLoaded + (testFirst-1);
camFrIndex = camFrIndex(camFrIndex >= camFirst);

% Ignore DIC images not requested for processing
dicFirst = dicLoaded + (testFirst-1);
dicFrIndex = 1:length(dicFiles);
dicFrIndex = dicFrIndex(dicFrIndex >= dicFirst)';

% Only keep images with both DIC and complete camera data
camtmp = camFrIndex-camFirst;
dictmp = dicFrIndex-dicFirst;
camFrIndex = (intersect(camtmp, dictmp)+camFirst);
dicFrIndex = (intersect(camtmp, dictmp)+dicFirst);

% Add the last DIC image without applied loading so that it can be
% subtracted from the images with applied loading (optional)
if ~isnan(SubIm)
    dicFrIndex = [dicLoaded+SubIm ; dicFrIndex];
end

dicData = dicDataLoad(dicFiles, dicFrIndex);

%% Proces results

% Number of rulers created during camera image processing
nRulers = max([max(camDataL(:,3)) max(camDataR(:,3))]);

VidSave = 0;

% Initiate waitbar
f = waitbar(1/length(dicFrIndex), 'Creating DIC images and video. Please wait...');

% Pre-define plot specifications
legendArray = {'\epsilon_{yy}' 'Disbond' 'Delamination_{1}' ...
    'Delamination_{2}' 'Delamination_{3}'};
lineSpec = {'--kx' '--ko' '--k^'};

if VidSave == 1
    % Create target folder for saving
    pathNew = fullfile(pathTarget, ['DIC [' num2str(filter) ' ' num2str(Ksize) ' ' num2str(CentPixExcl) ' ' num2str(SubIm)  '] ' datestr(now,30)]);
    mkdir(pathNew);
    
    % Create targetpath for video saving
    pathVideo = fullfile(pathNew, [clsName ' [' num2str(filter) ' ' num2str(Ksize) ' ' num2str(CentPixExcl) ' ' num2str(SubIm) '] ' datestr(now,30)]);
    
    % Prepare the new video file (4 frames/second)
    oVideo = VideoWriter(pathVideo, 'Uncompressed AVI');
    oVideo.FrameRate = 4;
    open(oVideo);
end

% Variable selection
var = dicData.ey;

for i = 60
    
    % Load cycle number
    nf = FirstImage+nStep*(dicFrIndex(i+(~isnan(SubIm)))-dicLoaded);
    
    % Extract data to be plotted
    if isnan(SubIm)
        rawZ = var(:,:,i);
        rawX = dicData.X(:,:,i);
        rawY = dicData.Y(:,:,i);
    elseif ~isnan(SubIm)
        rawZ = var(:,:,i+1)-var(:,:,1);
        rawX = dicData.X(:,:,i+1);
        rawY = dicData.Y(:,:,i+1);
    end
    
    % Adjust (X, Y) coordinates
    rawX_adj = rawX-(min(rawX(:))-rec(1))+dicEdgeGap;
    rawY_adj = rawY-(min(rawY(:))-rec(2))+dicEdgeGap;
    
    % Get the image frame data min and max values
    rawZmin = min(rawZ(~isnan(rawZ(:))));
    rawZmax = max(rawZ(~isnan(rawZ(:))));
    
    % Left and right disbond length
    bL = camDataL(camDataL(:,2) == camFrIndex(i),14);
    bR = camDataR(camDataR(:,2) == camFrIndex(i),14);
    if size(bL, 1) > size(bR, 1)
        % If more left camera images have been processed
        bLR = zeros(size(bL, 1), 2);
        bLR(:,1) = bL;
        bLR(1:size(bR, 1), 2) = bR;
    else
        % If more right camera images have been processed
        bLR = zeros(size(bR, 1),2);
        bLR(:, 2) = bR;
        bLR(1:size(bL, 1), 1) = bL;
    end
    
    % Interpolate rawZ over a mesh grid
    circX = circ(1)*cos(linspace(0, 2*pi, 200))+circ(2);
    circY = circ(1)*sin(linspace(0, 2*pi, 200))+circ(3);
    [intX, intY, intZ, intZi, intZii] = dicgriddata(rawX_adj, rawY_adj, rawZ, n(1), (2), 'fill', circX, circY);
    
    % Smooth the 1D results along y-axis (1-direction)
    intZs = nanconv(intZ, fspecial('average', 9), 'edge');
    intZis = smoothdata(intZi, 9, 'movmean', Ksize1D, 'omitnan');
    intZiis = smoothdata(intZii, 9, 'movmean', Ksize1D, 'omitnan');
    
    % Find the disbonded area
    [intZb, bA, intZp] = locateDisbondFront(intX, intY, intZs, 'findvalley', 'nofill', circX, circY);
    
    % Supress empty data points
    mask = (rawZ == 0);
    rawX_adj(mask) = NaN;
    rawY_adj(mask) = NaN;
    rawZ(mask) = NaN;
    
    % Create the DIC contour plot
    F = figure;
    title(['N_f = ' num2str(nf)]);
    pbaspect([2 6 1])
    
    subplot(1,3,1)
    hold on
    contourf(rawX_adj, rawY_adj, rawZ, linspace(rawZmin, rawZmax, 30), 'edgecolor', 'none');
    title('Raw data')
    colormap jet
    contourcbar
    % Include a rectangle to mark the edges of the overlap region
    rectangle('Position', rec, 'EdgeColor', 'black', 'LineWidth', 1.5);
    % Include the measured disbond/delamination positions
    for j = 1:size(bLR, 1)
        plot([rec(1) rec(1)+rec(3)], bLR(j, :)', lineSpec{j}, 'LineWidth', 2)
    end
    hold off
    caxis([rawZmin rawZmax])
    
    subplot(1,3,2)
    hold on
    contourf(intX, intY, intZ, 'edgecolor', 'none');
    title('Interpolated data')
    colormap jet
    contourcbar
    % Include a rectangle to mark the edges of the overlap region
    rectangle('Position', rec, 'EdgeColor', 'black', 'LineWidth', 1.5);
    % Include the measured disbond/delamination positions
    for j = 1:size(bLR, 1)
        plot([rec(1) rec(1)+rec(3)], bLR(j, :)', lineSpec{j}, 'LineWidth', 2)
    end
    hold off
    caxis([rawZmin rawZmax])
    
    subplot(1,3,3)
    hold on
    contourf(intX, intY, intZb, 'edgecolor', 'none');
    title('Disbonded area')
    colormap jet
    contourcbar
    % Include a rectangle to mark the edges of the overlap region
    rectangle('Position', rec, 'EdgeColor', 'black', 'LineWidth', 1.5);
    % Include the measured disbond/delamination positions
    for j = 1:size(bLR, 1)
        plot([rec(1) rec(1)+rec(3)], bLR(j, :)', lineSpec{j}, 'LineWidth', 2)
    end
    hold off
    caxis([rawZmin rawZmax])
    
    if VidSave == 1
        % Save figure to the new folder
        saveas(figure(i), fullfile(pathNew, [clsName '_Nf_' num2str(nf) '_[' num2str(filter) ' ' num2str(Ksize) ' ' num2str(CentPixExcl) ' ' num2str(SubIm) ']']), 'fig')
        saveas(figure(i), fullfile(pathNew, [clsName '_Nf_' num2str(nf) '_[' num2str(filter) ' ' num2str(Ksize) ' ' num2str(CentPixExcl) ' ' num2str(SubIm) ']']), 'png')
        saveas(figure(i), fullfile(pathNew, [clsName '_Nf_' num2str(nf) '_[' num2str(filter) ' ' num2str(Ksize) ' ' num2str(CentPixExcl) ' ' num2str(SubIm) ']']), 'eps')
        
        % Write each frame to the video file
        writeVideo(oVideo, getframe(figure(i)));
    end
    
    %close(F);
    
    waitbar(i/length(dicFrIndex), f);
end

if VidSave == 1
    close(oVideo);
end

% Close waitbar
close(f);
