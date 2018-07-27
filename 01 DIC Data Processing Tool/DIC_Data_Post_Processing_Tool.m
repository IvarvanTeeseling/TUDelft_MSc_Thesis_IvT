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

%% Collect and store data

% Select specimen
test_ID = {'00_01_Static_Reference';
    '00_02_Fatigue_Reference_30%';
    '00_03_Fatigue_Clamp_30%';
    '00_04_Fatigue_Reference_30%';
    '00_05_Fatigue_Bolt_30%'};
test_ID = test_ID{5};

% Get all .mat files
path_files  = fullfile(pwd, test_ID, '\00 DIC Data');
files       = dir(fullfile(path_files, '*.mat'));

% Read, extract and store all DIC data
for i = 1:length(files)
    % Load DIC data
    dic_rw(i) = load([files(i).folder '\' files(i).name]);
    
    % Pre-allocate memory
    if i == 1
        X   = zeros([size(dic_rw(i).X) length(files)]);
        Y   = zeros([size(dic_rw(i).X) length(files)]);
        W   = zeros([size(dic_rw(i).X) length(files)]);
        eyy = zeros([size(dic_rw(i).X) length(files)]);
        exx = zeros([size(dic_rw(i).x) length(files)]);
        e1  = zeros([size(dic_rw(i).X) length(files)]);
        e2  = zeros([size(dic_rw(i).X) length(files)]);
    end
    
    % Store data
    X(:,:,i)   = dic_rw(i).X;
    Y(:,:,i)   = dic_rw(i).Y;
    W(:,:,i)   = dic_rw(i).W;
    eyy(:,:,i) = dic_rw(i).eyy;
    exx(:,:,i) = dic_rw(i).exx;
    e1(:,:,i)  = dic_rw(i).e1;
    e2(:,:,i)  = dic_rw(i).e2;
end

% Suppress data points outside the area of interest
eyy(eyy==0) = NaN;
exx(exx==0) = NaN;
e1(e1==0)   = NaN;
e2(e2==0)   = NaN;
W(W==0)     = NaN;

% Temporary measured crack length
b = [25.4 25.4].*ones(length(files),2);

% Load cycles
N_start = 5000;
N_fig   = 9;
N_step  = 500;

%% Plot data - All DIC eyy contour plots & Video

Plot_1 = 1;

if Plot_1 == 1
    % New folder for saving
    path_nf = fullfile(pwd, test_ID, '01 Contour Plots', datestr(now,30));
    mkdir(path_nf);
    
    % New video path
    path_vd = fullfile(path_nf, test_ID);
    
    % Prepare the new video file (4 frames/second)
    vidObj = VideoWriter(path_vd);
    vidObj.FrameRate = 4;
    open(vidObj);
    
    % Variable to be plotted
    var = eyy;
    
    for i = 2:length(files)
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
        
        % Cycle number
        Nf = num2str(N_start+N_step*(i-N_fig));
        
        % Create plot
        figure(i)
        hold on
        contourf(X_tmp, Y_tmp, var_tmp, linspace(var_tmp_min, var_tmp_max, 20), 'edgecolor', 'none');
        colormap jet
        rectangle('Position', [-24 0 48 145], 'EdgeColor', 'black', 'LineWidth', 1);
        plot([-24 24], [b(i,1) b(i,2)], '--black', 'LineWidth', 1)
        pbaspect([1 6 1])
        contourcbar
        caxis([min(min(var(~isnan(var)))) max(max(var(~isnan(var))))])
        legend('\epsilon_{yy} [-]', 'Disbond position', 'location', 'northeastoutside')
        title(['N_f = ' Nf]);
        xlabel('z position [mm]')
        ylabel('y position [mm]')
        hold off
        
        % Save figure to the new folder
        saveas(figure(i), fullfile(path_nf, [test_ID '_Nf_' Nf]), 'fig')
        saveas(figure(i), fullfile(path_nf, [test_ID '_Nf_' Nf]), 'png')
        saveas(figure(i), fullfile(path_nf, [test_ID '_Nf_' Nf]), 'eps')
        
        % Write each frame to the video file
        writeVideo(vidObj, getframe(figure(i)));
        
        % Close the figure
        close(figure(i));
        
    end
    
    % Close the file.
    close(vidObj);
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
        title(['N_f = ' num2str(N_start+N_step*(i-N_fig)) ', \epsilon_{yy}']);
        xlabel('z position [mm]')
        ylabel('y position [mm]')
        hold off
    end
    
end
%% Plot data - Selected DIC eyy DIC contour vs. Analytical Model

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