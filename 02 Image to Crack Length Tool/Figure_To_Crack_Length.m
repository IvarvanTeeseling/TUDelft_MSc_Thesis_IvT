function varargout = Figure_To_Crack_Length(varargin)
% FIGURE_TO_CRACK_LENGTH MATLAB code for Figure_To_Crack_Length.fig
%      FIGURE_TO_CRACK_LENGTH, by itself, creates a new FIGURE_TO_CRACK_LENGTH or raises the existing
%      singleton*.
%
%      H = FIGURE_TO_CRACK_LENGTH returns the handle to a new FIGURE_TO_CRACK_LENGTH or the handle to
%      the existing singleton*.
%
%      FIGURE_TO_CRACK_LENGTH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGURE_TO_CRACK_LENGTH.M with the given input arguments.
%
%      FIGURE_TO_CRACK_LENGTH('Property','Value',...) creates a new FIGURE_TO_CRACK_LENGTH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Figure_To_Crack_Length_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Figure_To_Crack_Length_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Figure_To_Crack_Length

% Last Modified by GUIDE v2.5 13-Sep-2018 17:37:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Figure_To_Crack_Length_OpeningFcn, ...
    'gui_OutputFcn',  @Figure_To_Crack_Length_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Figure_To_Crack_Length is made visible.
function Figure_To_Crack_Length_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Figure_To_Crack_Length (see VARARGIN)

% Choose default command line output for Figure_To_Crack_Length
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Figure_To_Crack_Length wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Figure_To_Crack_Length_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbMakeMeasurement.
function pbMakeMeasurement_Callback(hObject, eventdata, handles)
% hObject    handle to pbMakeMeasurement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Fetch data
Ruler = getappdata(0, 'Ruler');
Files = getappdata(0, 'Files');
FileKey = getappdata(0, 'CurrentFileKey');
SF = getappdata(0, 'SF');
FJ = str2num(get(handles.etFigureJump, 'string'));
NpI = str2num(get(handles.etCyclesPerImage, 'string'));
NfI = str2num(get(handles.etFirstImageCycles, 'string'));

for i = 1:length(Ruler.hnd)
    % Get the ruler handle
    h = Ruler.hnd(i);
    
    % Check if still exists
    if Ruler.deleted(i) == 0
        % get the distance
        api = iptgetapi(h);
        dist = api.getDistance();
        pos = api.getPosition();
        
        if ~isappdata(0, 'Results')
            % Measurement counter
            cnt = 1;
            
            % Create matrix to store results
            Results = zeros(1,9);
            Results(1,1)    = 1;
            Results(1,2)    = FileKey;
            Results(1,3)    = Ruler.key(i);
            Results(1,4)    = pos(1,1);
            Results(1,5)    = pos(1,2);
            Results(1,6)    = pos(2,1);
            Results(1,7)    = pos(2,2);
            Results(1,8)    = dist;
            Results(1,9)    = SF.val;
            Results(1,10:14) = Results(1,4:8)*SF.val;
            Results(1,15)   = (FileKey-1)*NpI+NfI;
            
            % Store the results
            setappdata(0, 'Results', Results);
        else
            Results = getappdata(0, 'Results');
            
            % Update measurement counter
            if ~exist('cnt', 'var')
                cnt = max(Results(:,1))+1;
            end
            
            % Set results
            ResultsNew = [cnt ...
                FileKey ...
                Ruler.key(i)...
                pos(1,1) ...
                pos(1,2) ...
                pos(2,1) ...
                pos(2,2) ...
                dist ...
                SF.val ...
                pos(1,1)*SF.val ...
                pos(1,2)*SF.val ...
                pos(2,1)*SF.val ...
                pos(2,2)*SF.val ...
                dist*SF.val ...
                (FileKey-1)*NpI+NfI];
            
            % Overwrite existing file measurements
            if any(sum(Results(:,2:3)==[FileKey Ruler.key(i)], 2)==2)
                Results(sum(Results(:,2:3)==[FileKey Ruler.key(i)], 2)==2,:) = ResultsNew;
            else
                % Store the results
                Results = [Results ; ResultsNew];
            end
            % Set application data
            setappdata(0, 'Results', Results);
        end
        
    end
end

% Display results in the GUI
set(handles.tbResults, 'Data', Results);

% Only if multiple cells, not last figure and Auto next is toggled on
if iscell(Files.name) && FileKey+FJ <= length(Files.name) && get(handles.rbAutoNext, 'value') == 1
    % New figure whilst keeping the rulers and zoom settings
    set(Files.ImFig, 'CData', Files.read{FileKey+FJ})
    
    % Update File Key
    setappdata(0, 'CurrentFileKey', FileKey+FJ);
    
    % Update the GUI file overview
    set(handles.lbFileOverview, 'value', FileKey+FJ);
    set(handles.stCurrentFileKey, 'string', ['Curent file key: ' num2str(FileKey+FJ)]);
end

% --- Executes on button press in pbSelectImages.
function pbSelectImages_Callback(hObject, eventdata, handles)
% hObject    handle to pbSelectImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Ask for image files
[ImageFiles, FilePath] = uigetfile('Select the IMAGE FILE(s)','MultiSelect','on');

if isequal(ImageFiles,0) || isequal(FilePath,0)
    disp('User pressed cancel')
else
    % List file names in the GUI
    set(handles.lbFileOverview, 'string', ImageFiles);
    
    % Store file names and path
    Files.name = ImageFiles;
    Files.path = FilePath;
    
    % Initiate waitbar
    f = waitbar(1/length(Files.name),'Uploading files. Please wait...');
    
    % Read figures
    for i = 1:length(Files.name)
        % Read image file
        Files.read{i} = imread(fullfile(Files.path, Files.name{i}));
        
        % Update waitbar
        waitbar(i/length(Files.name),f);
    end
    
    % Close waitbar
    close(f);
    
    % Upload to application data
    setappdata(0, 'Files', Files);
end

% --- Executes on selection change in lbFileOverview.
function lbFileOverview_Callback(hObject, eventdata, handles)
% hObject    handle to lbFileOverview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lbFileOverview contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbFileOverview


% --- Executes during object creation, after setting all properties.
function lbFileOverview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbFileOverview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbStartMeasuring.
function pbStartMeasuring_Callback(hObject, eventdata, handles)
% hObject    handle to pbStartMeasuring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Files = getappdata(0, 'Files');

if ~iscell(Files.name)
    figure(1);
    Files.ImFig = imshow(Files.read);
else
    figure(1);
    Files.ImFig = imshow(Files.read{1});
end

set(handles.stCurrentFileKey, 'string', 'Curent file key: 1');

setappdata(0, 'Files', Files);
setappdata(0, 'CurrentFileKey', 1);


% --- Executes on button press in pbAddRuler.
function pbAddRuler_Callback(hObject, eventdata, handles)
% hObject    handle to pbAddRuler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make current figure for correct 'get current axes' (gca)
figure(1)

% Initiate a red ruler
RulerNew = imdistline(gca);

% Get application handler
api = iptgetapi(RulerNew);
api.setLabelVisible(false);

if isappdata(0, 'Ruler')
    % Collect current handles and keys
    Ruler = getappdata(0, 'Ruler');
    hnd = Ruler.hnd;
    key = Ruler.key;
    del = Ruler.deleted;
    
    % Add the new imdistline handle and key
    hnd = [hnd ; RulerNew];
    key = [key ; key(end)+1];
    del = [del ; 0];
    
    % Store
    Ruler.hnd = hnd;
    Ruler.key = key;
    Ruler.deleted = del;
    setappdata(0, 'Ruler', Ruler);
else
    % Create first handle, key and store
    Ruler.hnd = RulerNew;
    Ruler.key = 1;
    Ruler.deleted = 0;
    setappdata(0, 'Ruler', Ruler);
end

% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Ask for folder
folder = uigetdir();

% Create sepperate output for each ruler

Results = getappdata(0, 'Results');
Ruler   = getappdata(0, 'Ruler');

if ~isequal(folder,0)
    
    % Ask for file name
    prompt = {'Please give a file name (.xlsx will be added automatically):'};
    title = 'Input';
    answer = inputdlg(prompt,title);
    
    if isempty(answer)
        disp('Canceled saving.')
    else
        
        % Write to excel file in a new sheet
        xlswrite(fullfile(folder, answer{1}), Results, 'All Measurements', "A2");
        xlswrite(fullfile(folder, answer{1}), get(handles.tbResults, 'ColumnName')', 'All Measurements', "A1");
        
        for i = 1:size(Ruler.hnd, 1)
            % Find rows with corresponding ruler key
            ind = Results(:,3)==Ruler.key(i);
            
            % Isolate respective data
            res = Results(ind,:);
            
            xlswrite(fullfile(folder, answer{1}), res, ['Ruler ' num2str(Ruler.key(i))], "A2");
            xlswrite(fullfile(folder, answer{1}), get(handles.tbResults, 'ColumnName')', ['Ruler ' num2str(Ruler.key(i))], "A1");
        end
    end
end


% --- Executes on button press in pbSetScaling.
function pbSetScaling_Callback(hObject, eventdata, handles)
% hObject    handle to pbSetScaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer = questdlg('After pressing CONTINUE, please position the red ruler and press ANY KEY. You will be asked for the known dimensions next.', ...
    'Instructions', ...
    'Continue','Cancel','Cancel');

if strcmp(answer, 'Continue')
    % Make current figure for correct 'get current axes' (gca)
    figure(1)
    
    % Initiate a red ruler
    h = imdistline(gca);
    setColor(h, 'red');
    api = iptgetapi(h);
    api.setLabelVisible(false);
    
    % Wait for any key
    pause()
    
    % Ask for known dimensions
    prompt = {'Please position the red ruler correctly and enter the known dimension:','Enter unit:'};
    title = 'Input';
    dims = [1 35];
    definput = {'1','mm'};
    answer = inputdlg(prompt,title,dims,definput);
    
    if ~isempty(answer)
        % Extract the known dimension
        KnownDim = str2double(answer{1});
        
        % Ruler distance and position
        SF.dist = api.getDistance();
        SF.pos = api.getPosition();
        
        % Scale Factor
        SF.val = KnownDim/SF.dist;
        
        % Unit
        SF.uni = answer{2};
        
        % Store results
        setappdata(0, 'SF', SF);
        
        % Update static text
        str = [num2str(SF.val) ' [' answer{2} '] '];
        set(handles.stSF, 'string', str);
    end
    
    % Delete the ruler
    delete(h);
    
end

% --- Executes during object creation, after setting all properties.
function tbResults_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tbResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Set column names
set(hObject, 'columnname', {'Measure Index', 'File Index', 'Ruler Index', 'x_0', 'x_1', 'y_0', 'y_1', 'l', 'SF', 'X_0', 'X_1', 'Y_0', 'Y_1', 'L', 'Load Cycle'});


% --- Executes on button press in pbDeleteRuler.
function pbDeleteRuler_Callback(hObject, eventdata, handles)
% hObject    handle to pbDeleteRuler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isappdata(0, 'Ruler')
    % Get rulers
    Ruler = getappdata(0, 'Ruler');
    
    % Ask for ruler to be deleted
    prompt = {'Enter ruler to be deleted:'};
    title = 'Input';
    dims = [1 35];
    definput = {'0'};
    answer = str2double(inputdlg(prompt,title,dims,definput));
    
    % Continue of ruler exists
    if answer <= length(Ruler.hnd) && answer > 0 && Ruler.deleted(answer) == 0
        delete(Ruler.hnd(answer));
        Ruler.deleted(answer) = 1;
        
        setappdata(0, 'Ruler', Ruler)
    end
end

% --- Executes during object creation, after setting all properties.
function pbClear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pbClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pbNextImages.
function pbNextImages_Callback(hObject, eventdata, handles)
% hObject    handle to pbNextImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Files = getappdata(0, 'Files');
FileKey = getappdata(0, 'CurrentFileKey');

% Only if multiple figures and not the last one
if iscell(Files.name) && FileKey < length(Files.name)
    % New figure whilst keeping the rulers and zoom settings
    set(Files.ImFig, 'CData', Files.read{FileKey+1})
    
    % Update File Key
    setappdata(0, 'CurrentFileKey', FileKey+1);
    
    % Update the GUI file overview
    set(handles.lbFileOverview, 'value', FileKey+1);
    set(handles.stCurrentFileKey, 'string', ['Curent file key: ' num2str(FileKey+1)]);
end

% --- Executes on button press in pbPreviousImage.
function pbPreviousImage_Callback(hObject, eventdata, handles)
% hObject    handle to pbPreviousImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Files = getappdata(0, 'Files');
FileKey = getappdata(0, 'CurrentFileKey');

% Only if multiple figures and not the last one
if iscell(Files.name) && FileKey > 1
    % New figure whilst keeping the rulers and zoom settings
    set(Files.ImFig, 'CData', Files.read{FileKey-1})
    
    % Update File Key
    setappdata(0, 'CurrentFileKey', FileKey-1);
    
    % Update the GUI file overview
    set(handles.lbFileOverview, 'value', FileKey-1);
    set(handles.stCurrentFileKey, 'string', ['Curent file key: ' num2str(FileKey-1)]);
end

% --- Executes on button press in rbAutoNext.
function rbAutoNext_Callback(hObject, eventdata, handles)
% hObject    handle to rbAutoNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbAutoNext



function etImageJumping_Callback(hObject, eventdata, handles)
% hObject    handle to etImageJumping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etImageJumping as text
%        str2double(get(hObject,'String')) returns contents of etImageJumping as a double


% --- Executes during object creation, after setting all properties.
function etImageJumping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etImageJumping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etCyclesPerImage_Callback(hObject, eventdata, handles)
% hObject    handle to etCyclesPerImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etCyclesPerImage as text
%        str2double(get(hObject,'String')) returns contents of etCyclesPerImage as a double


% --- Executes during object creation, after setting all properties.
function etCyclesPerImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etCyclesPerImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etFirstImageCycles_Callback(hObject, eventdata, handles)
% hObject    handle to etFirstImageCycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etFirstImageCycles as text
%        str2double(get(hObject,'String')) returns contents of etFirstImageCycles as a double


% --- Executes during object creation, after setting all properties.
function etFirstImageCycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etFirstImageCycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function etFigureJump_Callback(hObject, eventdata, handles)
% hObject    handle to etFigureJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etFigureJump as text
%        str2double(get(hObject,'String')) returns contents of etFigureJump as a double


% --- Executes during object creation, after setting all properties.
function etFigureJump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etFigureJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pbClear.
function pbClear_Callback(hObject, eventdata, handles)
% hObject    handle to pbClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Ask for confirmation
promptMessage = sprintf('Do you want to Continue? All data will be lost!');
button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');

% Only delete when confirmed
if strcmpi(button, 'Continue')
    if isappdata(0, 'Ruler')
        rmappdata(0, 'Ruler')
    end
    
    if isappdata(0, 'Ruler')
        rmappdata(0, 'Ruler')
    end
    if isappdata(0, 'SF')
        rmappdata(0, 'SF')
    end
    
    if isappdata(0, 'CurrentFileKey')
        rmappdata(0, 'CurrentFileKey')
    end
    
    if isappdata(0, 'Results')
        rmappdata(0, 'Results')
    end
    
    % Clear results overview
    set(handles.tbResults, 'Data', []);
    
    % Close the figure
    close(figure(1));
end

% --- Executes on button press in pbCloseFigure.
function pbCloseFigure_Callback(hObject, eventdata, handles)
% hObject    handle to pbCloseFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isappdata(0, 'Ruler')
    rmappdata(0, 'Ruler')
end

if isappdata(0, 'Files')
    rmappdata(0, 'Files')
end

if isappdata(0, 'Ruler')
    rmappdata(0, 'Ruler')
end

if isappdata(0, 'SF')
    rmappdata(0, 'SF')
end

if isappdata(0, 'CurrentFileKey')
    rmappdata(0, 'CurrentFileKey')
end

if isappdata(0, 'Results')
    rmappdata(0, 'Results')
end

% Close the figure
close('Figure_To_Crack_Length')


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)  


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
