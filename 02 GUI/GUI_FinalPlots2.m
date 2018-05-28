function varargout = GUI_FinalPlots2(varargin)
% GUI_FINALPLOTS2 MATLAB code for GUI_FinalPlots2.fig
%      GUI_FINALPLOTS2, by itself, creates a new GUI_FINALPLOTS2 or raises the existing
%      singleton*.
%
%      H = GUI_FINALPLOTS2 returns the handle to a new GUI_FINALPLOTS2 or the handle to
%      the existing singleton*.
%
%      GUI_FINALPLOTS2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FINALPLOTS2.M with the given input arguments.
%
%      GUI_FINALPLOTS2('Property','Value',...) creates a new GUI_FINALPLOTS2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FinalPlots2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FinalPlots2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FinalPlots2

% Last Modified by GUIDE v2.5 25-May-2018 15:52:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FinalPlots2_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FinalPlots2_OutputFcn, ...
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


% --- Executes just before GUI_FinalPlots2 is made visible.
function GUI_FinalPlots2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FinalPlots2 (see VARARGIN)

% Choose default command line output for GUI_FinalPlots2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_FinalPlots2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FinalPlots2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_Cycle_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Cycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

x           = getappdata(0, 'x');
MinorSum    = getappdata(0, 'MinorSum');
N           = getappdata(0, 'N');
dbdN        = getappdata(0, 'dbdN');

Sa_nom_ad1  = getappdata(0, 'Sa_nom_ad1');
Sm_nom_ad1  = getappdata(0, 'Sm_nom_ad1');
St_nom_ad1  = Sa_nom_ad1 + Sm_nom_ad1;

value = get(hObject, 'value');
value = round(value);

% Update static text
set(handles.text_b,'String',num2str(x(1,value)))
set(handles.text_N,'String',num2str(N(value)*1e-3))

% Update the plots
set(handles.PlotMinor1,'XData',x(1,:),'YData',MinorSum(value,:));
set(handles.PlotMinor2,'XData',x(1,:),'YData',ones(1,size(x,2)));

set(handles.PlotdbdN,'XData',N(1:value),'YData',dbdN(1:value));

set(handles.PlotSm_ad1,'XData',x(value, value:end),'YData',Sm_nom_ad1(value, value:end,1));
set(handles.PlotStot_ad1,'XData',x(value, value:end),'YData',St_nom_ad1(value, value:end,1));

x = x(1,1:size(x,1))-x(1,1);
set(handles.Plotb,'XData',N(1:value),'YData',x(1,1:value));


% --- Executes during object creation, after setting all properties.
function slider_Cycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Cycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

x = getappdata(0, 'x');

set(hObject, 'Min', 1);
set(hObject, 'Max', size(x,1));
set(hObject, 'Value', size(x,1));
set(hObject, 'SliderStep', [10/(size(x,1)-1) , 100/(size(x,1)-1)]);


% --- Executes during object creation, after setting all properties.
function axes_Minor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_Minor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_Minor

x = getappdata(0, 'x');
MinorSum = getappdata(0, 'MinorSum');

hold on
handles.PlotMinor1 = plot(x(1,:), MinorSum(end,:),'b');
handles.PlotMinor2 = plot(x(1,:), ones(1,size(x,2)),'r');
hold off
title('Fatigue Accumulation')
ylabel('Minor Damage Index [-]')
xlabel('Distance x [mm]')
xlim([x(1,1) x(1,end)]);
if max(MinorSum(:)) > 1
    ylim([0 max(MinorSum(:))]);
else
    ylim([0 1]);
end
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_b
% if isappdata(0, 'n')
%     n = getappdata(0, 'n');
% else
%     n = evalin('base','n');
%     setappdata(0, 'n', n);
% end
% 
% if isappdata(0, 'x')
%     x = getappdata(0, 'x');
% else
%     x = evalin('base','x00');
%     setappdata(0, 'x', x);
% end

N = getappdata(0, 'N');
x = getappdata(0, 'x');

x = x(1,1:size(x,1))-x(1,1);

hold on
handles.Plotb = plot(N, x,'b');
hold off
title('Crack length (b)')
ylabel('Length [mm]')
xlabel('Number of load cycles (N) [-]')
xlim([0 N(end)]);
ylim([0 x(end)]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_dbdN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_dbdN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_dbdN

% if isappdata(0, 'n')
%     n = getappdata(0, 'n');
% else
%     n = evalin('base','n');
%     setappdata(0, 'n', n);
% end
% 
% if isappdata(0, 'dbdN')
%     dbdN = getappdata(0, 'dbdN');
% else
%     dbdN = evalin('base','dbdN');
%     setappdata(0, 'dbdN', dbdN);
% end

N = getappdata(0, 'N');
dbdN = getappdata(0, 'dbdN');

hold on
handles.PlotdbdN = plot(N, dbdN,'b');
hold off
title('Crack growth rate (db/dN)')
ylabel('Growth Rate [mm/cycle]')
xlabel('Number of load cycles (N) [-]')
xlim([0 N(end)]);
ylim([0 max(dbdN)]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes_StressCycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_StressCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_StressCycle

x = getappdata(0, 'x');
Sa_nom_ad1 = getappdata(0, 'Sa_nom_ad1');
Sm_nom_ad1 = getappdata(0, 'Sm_nom_ad1');
Sy         = getappdata(0, 'Sy');
St_nom_ad1 = Sa_nom_ad1 + Sm_nom_ad1;

hold on
handles.PlotSm_ad1   = plot(x(1,:), Sm_nom_ad1(end,:),'r');
handles.PlotStot_ad1 = plot(x(1,:), St_nom_ad1(end,:),'b');
handles.PlotSy       = plot(x(1,:), Sy*ones(1,size(x,2)),'g');
hold off
title('Load Cycle Stress')
legend('Mean S_{m,nom}', 'Total Stress S_{tot}','Yield Stress \sigma_{y}')
ylabel('Stress [MPa]')
xlabel('Distance x [mm]')
xlim([x(1,1) x(1,end)]);
if max(St_nom_ad1(:)) >= Sy
    ylim([min(Sm_nom_ad1(:)) max(St_nom_ad1(:))*1.1]);
else
    ylim([min(Sm_nom_ad1(:)) Sy*1.1]);
end
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_G_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_G (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_G

x   = getappdata(0, 'x');
GI  = getappdata(0, 'GI');
GII = getappdata(0, 'GII');
G   = GI + GII;

hold on
handles.PlotG   = plot(x(1,1:size(x,1)), G(:,:,2),'r');
handles.PlotGI  = plot(x(1,1:size(x,1)), GI(:,:,2),'b');
handles.PlotGII = plot(x(1,1:size(x,1)), GII(:,:,2),'g');
hold off
title('Strain Energy Release Rate')
legend('G', 'G_{I}','G_{II}')
ylabel('SERR')
xlabel('Distance x [mm]')
xlim([x(1,1) x(1,end)]);
ylim([0 max(G(:))]);
grid on

% Update handles structure
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function axes_Geq_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes_Geq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes_Geq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_Geq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_Geq

x   = getappdata(0, 'x');
dG1_eq  = getappdata(0, 'dG1_eq');

hold on
handles.PlotdG1_eq = plot(x(1,1:size(x,1)), dG1_eq (:,:),'r');
hold off
title('Equivalent SERR Range')
legend('dG_{1,eq}')
ylabel('dSERR')
xlabel('Distance x [mm]')
xlim([x(1,1) x(1,end)]);
ylim([0 max(dG1_eq (:))]);
grid on

% Update handles structure
guidata(hObject, handles);
