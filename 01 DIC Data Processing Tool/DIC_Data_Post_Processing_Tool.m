%% Description
%
%
%
%
%
%
%

%% Initialize MATLAB
clear; close all; clc

% Set timer
tic

%% File selection - 01_03_Fat_Ref_30kN > No taper

% Camera files
path    = 'D:\Thesis Experimental Data\Run_01_CLS\01_03_Fat_Ref_30kN\Data Processing';
% Note: right and left seemed to be switched for this specimen
nameL   = 'Camera Post Processing - 01_03_Fat_Ref_30kN (Right Camera).txt';
nameR   = 'Camera Post Processing - 01_03_Fat_Ref_30kN (Left Camera).txt';

% DIC files
path_files  = fullfile('D:\Thesis Experimental Data\Run_01_CLS\01_03_Fat_Ref_30kN\DIC\MATLAB Data');
path_store  = fullfile('D:\Thesis Experimental Data\Run_01_CLS\01_03_Fat_Ref_30kN\Data Processing');
SpecName    = '01_03_Fat_Ref_30kN';

% Load cycles
N_start     = 500;
N_step      = 500;

% DIC Image selection
FInitialDIC  = 7;

%% File selection - 01-07-Fat_Ref_30kN > Stepped taper

% % Camera files
% path    = 'D:\Thesis Experimental Data\Run_01_CLS\01_07_Fat_Ref_30kN\Data Processing';
% % Note: no right camera for this specimen; was pointed at the tapered edge
% % so monitor Al fatige crack initiation
% nameL   = 'Crack Length - Left Camera - 01_07_Fat_Ref_30kN.txt';
% nameR   = nameL;
% 
% % DIC files
% path_files  = fullfile('D:\Thesis Experimental Data\Run_01_CLS\01_07_Fat_Ref_30kN\DIC\MATLAB Data');
% path_store  = fullfile('D:\Thesis Experimental Data\Run_01_CLS\01_07_Fat_Ref_30kN\Data Processing');
% SpecName    = '01_07_Fat_Ref_30kN';
% 
% % Load cycles
% N_start     = 500;
% N_step      = 500;
% 
% % DIC Image selection
% FInitialDIC  = 1;

%% Processing

% Get all files
files   = dir(fullfile(path_files, '*.mat'));
bLeft   = textread(fullfile(path, nameL));
bRight  = textread(fullfile(path, nameR));

% Frames analyzed
if max(bLeft(:,2)) > max(bRight(:,2))
    FrameIndex = unique(bRight(:,2));
else
    FrameIndex = unique(bLeft(:,2));
end

% Number of rulers
nRuler = max([max(bLeft(:,3)) max(bRight(:,3))]);

% Analyzed frames
nFrame = length(FrameIndex);

% Initiate waitbar
f = waitbar(1/nFrame,'Uploading MATLAB data. Please wait...');

% Read, extract and store all DIC data
for i = 1:nFrame
    
    ind = FrameIndex(i);
    
    % Load DIC data
    dic_rw(i) = load(fullfile(files(ind+FInitialDIC-1).folder, files(ind+FInitialDIC-1).name));
    
    % Pre-allocate memory
    if i == 1
        X   = zeros([size(dic_rw(i).X) nFrame]);
        Y   = zeros([size(dic_rw(i).X) nFrame]);
        W   = zeros([size(dic_rw(i).X) nFrame]);
        eyy = zeros([size(dic_rw(i).X) nFrame]);
        exx = zeros([size(dic_rw(i).x) nFrame]);
        e1  = zeros([size(dic_rw(i).X) nFrame]);
        e2  = zeros([size(dic_rw(i).X) nFrame]);
    end
    
    % Store data
    X(:,:,i)   = dic_rw(i).X;
    Y(:,:,i)   = dic_rw(i).Y;
    W(:,:,i)   = dic_rw(i).W;
    eyy(:,:,i) = dic_rw(i).eyy;
    exx(:,:,i) = dic_rw(i).exx;
    e1(:,:,i)  = dic_rw(i).e1;
    e2(:,:,i)  = dic_rw(i).e2;
    
    % Update waitbar
    waitbar(i/nFrame,f);
end

% Close waitbar
close(f);

% Suppress data points outside the area of interest
eyy(eyy==0) = NaN;
exx(exx==0) = NaN;
e1(e1==0)   = NaN;
e2(e2==0)   = NaN;
W(W==0)     = NaN;

%% Plot data - All DIC eyy contour plots & Video

Plot_1 = 1;

% Legend array
LArray = {'\epsilon_{yy}' 'Disbond' 'Delamination_{1}' 'Delamination_{2}'};

% Line spec
LSpec = {'--kx' '--ko' '--k^'};

if Plot_1 == 1
    
    % Initiate waitbar
    f = waitbar(1/nFrame,'Creating DIC images and video. Please wait...');
    
    % New folder for saving
    path_nf = fullfile(path_store, ['DIC Processing ' datestr(now,30)]);
    mkdir(path_nf);
    
    % New video path
    path_vd = fullfile(path_nf, SpecName);
    
    % Prepare the new video file (4 frames/second)
    vidObj = VideoWriter(path_vd, 'Uncompressed AVI');
    vidObj.FrameRate = 4;
    open(vidObj);
    
    % Variable to be plotted
    var = eyy;
    
    for i = 1:nFrame
        % Extract data to be plotted
        var_tmp     = var(:,:,i);
        X_tmp       = X(:,:,i);
        Y_tmp       = Y(:,:,i);
        
        % Suppress empty data points
        X_tmp(isnan(var_tmp)) = NaN;
        Y_tmp(isnan(var_tmp)) = NaN;
        
        % Data set minimum and maximum values
        var_tmp_min = min(var_tmp(~isnan(var_tmp(:))));
        var_tmp_max = max(var_tmp(~isnan(var_tmp(:))));
        
        % Adjust (x, y) coordinates
        Y_tmp = Y_tmp - (max(Y_tmp(:)));
        
        % Cycle number
        Nf = N_start+N_step*(FrameIndex(i)-1);
        
        % Get measured crack positions
        bL = bLeft(bLeft(:,2)==FrameIndex(i), 14);
        bR = bRight(bRight(:,2)==FrameIndex(i), 14);
        if size(bL, 1) > size(bR, 1)
            bLR                 = zeros(size(bL,1),2);
            bLR(:,1)            = -1*bL;
            bLR(1:size(bR,1),2) = -1*bR;
        else
            bLR                 = zeros(size(bR,1),2);
            bLR(:,2)            = -1*bR;
            bLR(1:size(bL,1),1) = -1*bL;
        end
        
        % Create plot
        figure(i)
        hold on
        contourf(X_tmp, Y_tmp, var_tmp, linspace(var_tmp_min, var_tmp_max, 30), 'edgecolor', 'none');
        colormap jet
        rectangle('Position', [-22.2 -145 44.4 145], 'EdgeColor', 'black', 'LineWidth', 1.5);
        for j = 1:size(bLR,1)
            plot([-22.2 22.2], bLR(j,:)', LSpec{j}, 'LineWidth', 2)
        end
        pbaspect([2 6 1])
        contourcbar
        caxis([var_tmp_min var_tmp_max])
        legend(LArray{1:size(bLR,1)+1}, 'location', 'northeastoutside')
        title(['N_f = ' num2str(Nf)]);
        xlabel('z position [mm]')
        ylabel('y position [mm]')
        hold off
        
        % Save figure to the new folder
        saveas(figure(i), fullfile(path_nf, [SpecName '_Nf_' num2str(Nf)]), 'fig')
        saveas(figure(i), fullfile(path_nf, [SpecName '_Nf_' num2str(Nf)]), 'png')
        saveas(figure(i), fullfile(path_nf, [SpecName '_Nf_' num2str(Nf)]), 'eps')
        
        % Write each frame to the video file
        writeVideo(vidObj, getframe(figure(i)));
        
        % Close the figure
        close(figure(i));
        
        % Update waitbar
        waitbar(i/nFrame,f);
        
    end
    
    % Close the file.
    close(vidObj);
    
    % Close waitbar
    close(f);
end

%% Plot data - Subplot eyy countour for selected DIC foto's

Plot_2 = 0;

if Plot_2 == 1
    % Variable to be plotted
    var     = eyy;
    
    % Select foto's
    select_vec = [2 5 10 15 30 45 150];
    for j = 1:length(select_vec)
        % DIC data index
        i = select_vec(j);
        
        % Extract data to be plotted
        var_tmp     = var(:,:,i)-var(:,:,2);
        X_tmp       = X(:,:,i);
        Y_tmp       = Y(:,:,i)+abs(min(min(Y(:,:,i))));
        
        % Suppress empty data points
        X_tmp(isnan(var_tmp)) = NaN';
        Y_tmp(isnan(var_tmp)) = NaN';
        
        % Data set minimum and maximum values
        var_tmp_min = min(var_tmp(~isnan(var_tmp(:))));
        var_tmp_max = max(var_tmp(~isnan(var_tmp(:))));
        
        % Create plot
        figure(length(select_vec)+1)
        subplot(1,length(select_vec),j)
        hold on
        contourf(X_tmp, Y_tmp, var_tmp, linspace(var_tmp_min,var_tmp_max,20),'edgecolor','none');
        rectangle('Position',[-24 0 48 145],'EdgeColor','black', 'LineWidth', 1);
        plot([-24 24], [b(i,1) b(i,2)], ':black', 'LineWidth', 2)
        contourcbar
        caxis([min(min(var(~isnan(var)))) max(max(var(~isnan(var))))])
        colormap jet
        title(['N_f = ' num2str(N_start+N_step*(i-FInitialDIC)) ', \epsilon_{yy}']);
        xlabel('z position [mm]')
        ylabel('y position [mm]')
        hold off
    end
    
end
%% Plot data - Selected DIC eyy DIC contour vs. Analytical Model

% n = 12;
% z = eyy(:,:,n);
% x = X(:,:,n);
% y = Y(:,:,n);
% 
% z(isnan(z)) = 0;
% x(isnan(x)) = 0;
% y(isnan(y)) = 0;
% 
% % x(isnan(z)) = NaN;
% yy = y;
% yy(z==0) = NaN;
% y = y - (max(yy(:)));
% 
% yq = linspace(0, -120);
% xq = zeros(size(yq));
% 
% z_int = griddata(x, y, z, xq, yq);


Plot_3 = 0;

if Plot_3 == 1
    
    load('variables.mat');
    
    % Extract data to be plotted
    var_tmp     = eyy(:,:,3);
    X_tmp       = X(:,:,3);
    Y_tmp       = Y(:,:,3)+abs(min(min(Y(:,:,3))));
    
    % Suppress empty data points
    X_tmp(isnan(var_tmp)) = NaN';
    Y_tmp(isnan(var_tmp)) = NaN';
    
    % Data set minimum and maximum values
    var_tmp_min = min(var_tmp(~isnan(var_tmp(:))));
    var_tmp_max = max(var_tmp(~isnan(var_tmp(:))));
    
    z_int = griddata(x,(y+0.0254)*1000,z,X_tmp,Y_tmp);
    
    if min(z(:)) < var_tmp_min
        var_tmp_min = min(z(:));
    end
    
    if max(z(:)) > var_tmp_max
        var_tmp_max = max(z(:));
    end
    
    % Remove singularities
    err = (z_int-var_tmp)./var_tmp*100;
    err(abs(err)>200) = NaN;
    
    figure(length(select_vec)+2)
    subplot(1,3,1)
    hold on
    contourf(x, (y+0.0254)*1000, z, linspace(min(z(:)),max(z(:)),25),'edgecolor','none');
    rectangle('Position',[-24 0 48 145])
    contourcbar
    title('\epsilon_{yy} - Model prediction')
    xlabel('x location')
    ylabel('y location')
    caxis([var_tmp_min var_tmp_max]);
    hold off
    subplot(1,3,2)
    hold on
    contourf(X_tmp, Y_tmp, var_tmp, linspace(var_tmp_min,var_tmp_max,25),'edgecolor','none');
    rectangle('Position',[-24 0 48 145])
    contourcbar
    title('\epsilon_{yy} - Experimental data')
    xlabel('x location')
    ylabel('y location')
    caxis([var_tmp_min var_tmp_max]);
    hold off
    % subplot(1,3,3)
    % hold on
    % contourf(X_tmp, Y_tmp,err, linspace(min(err(:)),max(err(:)),25),'edgecolor','none');
    % rectangle('Position',[-24 0 48 145])
    % contourcbar
    % title('\epsilon_{yy} - Experimental data')
    % xlabel('x location')
    % ylabel('y location')
    % hold off
    
end

% End timer
toc