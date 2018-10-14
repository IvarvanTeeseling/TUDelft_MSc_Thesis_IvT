function varargout = GUI_FinalPlots3(varargin)
% GUI_FINALPLOTS3 MATLAB code for GUI_FinalPlots3.fig
%      GUI_FINALPLOTS3, by itself, creates a new GUI_FINALPLOTS3 or raises the existing
%      singleton*.
%
%      H = GUI_FINALPLOTS3 returns the handle to a new GUI_FINALPLOTS3 or the handle to
%      the existing singleton*.
%
%      GUI_FINALPLOTS3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FINALPLOTS3.M with the given input arguments.
%
%      GUI_FINALPLOTS3('Property','Value',...) creates a new GUI_FINALPLOTS3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FinalPlots3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FinalPlots3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FinalPlots3

% Last Modified by GUIDE v2.5 13-Aug-2018 17:11:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FinalPlots3_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FinalPlots3_OutputFcn, ...
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


% --- Executes just before GUI_FinalPlots3 is made visible.
function GUI_FinalPlots3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FinalPlots3 (see VARARGIN)

% Choose default command line output for GUI_FinalPlots3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_FinalPlots3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FinalPlots3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

xAB         = getappdata(0, 'xAB');
Sa_nom_ACB  = getappdata(0, 'Sa_nom_ACB');
Sm_nom_ACB  = getappdata(0, 'Sm_nom_ACB');
Sy          = getappdata(0, 'Sy');
Su          = getappdata(0, 'Su');
St_nom_ACB  = Sa_nom_ACB+Sm_nom_ACB;
b           = getappdata(0, 'b');
l_A         = getappdata(0, 'l_A');
b           = b+l_A(1);

hold on
handles.PlotSt_ACBt     = plot(xAB(1,:), St_nom_ACB(end,:,2), 'r');
handles.PlotSt_ACBl     = plot(xAB(1,:), St_nom_ACB(end,:,1), 'b');
handles.PlotSy          = plot(xAB(1,:), Sy*ones(1,size(xAB,2)), 'g');
handles.PlotSu          = plot(xAB(1,:), Su*ones(1,size(xAB,2)), 'c');
xlim([xAB(1,1) xAB(1,end)]);
ylim([min(Sm_nom_ACB(:)) Su*1.5]);
handles.PlotDisbondPos3 = plot([b(end) b(end)], ylim, '-mo');
hold off
grid on
title('Total Ply Stresses')
legend('S_{T,xx,max}', 'S_{L,xx,max}', '\sigma_{y}' , '\sigma_{u}', 'Disbond [mm]', 'location', 'southwest')
ylabel('Stress [MPa]')
xlabel('x_{AB} [mm]')

% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

xAB         = getappdata(0, 'xAB');
Minor       = getappdata(0, 'Minor');
dMinor      = getappdata(0, 'dMinor');
N           = getappdata(0, 'N');
dbdN        = getappdata(0, 'dbdN');
b           = getappdata(0, 'b');
l_A         = getappdata(0, 'l_A');

Sa_nom_ACB  = getappdata(0, 'Sa_nom_ACB');
Sm_nom_ACB  = getappdata(0, 'Sm_nom_ACB');
St_nom_ACB  = Sa_nom_ACB + Sm_nom_ACB;

value = get(hObject, 'value');
value = round(value);

% Update static text
set(handles.stCrackLength, 'String', num2str(b(value)));
set(handles.stLoadCycle, 'String', num2str(N(value)));

% Update the plots
set(handles.PlotMinort, 'XData', xAB(1,:), 'YData', Minor(value,:,2));
set(handles.PlotMinorl, 'XData', xAB(1,:), 'YData', Minor(value,:,1));

set(handles.PlotdMinort, 'XData', xAB(1,:), 'YData', dMinor(value,:,2));
set(handles.PlotdMinorl, 'XData', xAB(1,:), 'YData', dMinor(value,:,1));

set(handles.PlotdbdN, 'XData', N(1:value), 'YData', dbdN(1:value));

set(handles.PlotSt_ACBt, 'XData', xAB, 'YData', St_nom_ACB(value, :, 2));
set(handles.PlotSt_ACBl, 'XData', xAB, 'YData', St_nom_ACB(value, :, 1));

set(handles.PlotDisbondPos1, 'XData', [b(value)+l_A(1) b(value)+l_A(1)]);
set(handles.PlotDisbondPos2, 'XData', [b(value)+l_A(1) b(value)+l_A(1)]);
set(handles.PlotDisbondPos3, 'XData', [b(value)+l_A(1) b(value)+l_A(1)]);

set(handles.Plotb,'XData', N(1:value), 'YData', b(1:value));

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

b = getappdata(0, 'b');

set(hObject, 'Min', 1);
set(hObject, 'Max', length(b));
set(hObject, 'Value', length(b));
set(hObject, 'SliderStep', [10/(length(b)-1) , 100/(length(b)-1)]);

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2

N = getappdata(0, 'N');
b = getappdata(0, 'b');

hold on
handles.Plotb = plot(N(1:length(b)), b, 'b');
hold off
title('Crack length (b)')
ylabel('Length [mm]')
xlabel('N [-]')
xlim([0 N(end)]);
ylim([0 b(end)]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes5

x       = getappdata(0, 'x');
GI      = getappdata(0, 'GI');
GII     = getappdata(0, 'GII');
G       = GI + GII;
dG_eq   = getappdata(0, 'dG1_eq');

hold on
handles.PlotG       = plot(x(1,1:size(x,1)), G(:,:,2),'r');
handles.PlotGI      = plot(x(1,1:size(x,1)), GI(:,:,2),'b');
handles.PlotGII     = plot(x(1,1:size(x,1)), GII(:,:,2),'g');
handles.PlotdG_eq   = plot(x(1,1:size(x,1)), dG_eq(:,:),'c');
hold off
title('Strain Energy Release Rate')
legend('G', 'G_{I}','G_{II}', 'dG_{I,eq}')
ylabel('SERR [J/m^2]')
xlabel('x_{BB} [mm]')
xlim([x(1,1) x(1,end)]);
ylim([0 max([max(G(:)) max(dG_eq(:))])]);
grid on

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4

N       = getappdata(0, 'N');
dbdN    = getappdata(0, 'dbdN');

hold on
handles.PlotdbdN = plot(N, dbdN,'b');
hold off
title('Crack growth rate (db/dN)')
ylabel('Growth Rate [mm/cycle]')
xlabel('N [-]')
xlim([0 N(end)]);
ylim([0 max(dbdN)]);
grid on

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3

Minor   = getappdata(0, 'Minor');
xAB     = getappdata(0, 'xAB');
b       = getappdata(0, 'b');
l_A     = getappdata(0, 'l_A');
b       = b+l_A(1);

hold on
handles.PlotMinort  = plot(xAB(1,:), Minor(end, :, 2), 'r');
handles.PlotMinorl  = plot(xAB(1,:), Minor(end, :, 1), 'b');
plot([0 xAB(1,end)], [1 1], 'g')
xlim([0 xAB(1,end)]);
ylim([0 max(Minor(:))]);
handles.PlotDisbondPos2 = plot([b(end) b(end)], ylim, '-mo');
hold off
title('Minor index')
xlabel('x_{AB} [mm]')
ylabel('Minor index [-]')
legend('Top Ply', 'Bottom ply', 'Initiation threshold', 'Disbond [mm]')
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3

dMinor  = getappdata(0, 'dMinor');
xAB     = getappdata(0, 'xAB');
b       = getappdata(0, 'b'); 
l_A     = getappdata(0, 'l_A');
b       = b+l_A(1);

hold on
handles.PlotdMinort = plot(xAB(1,:), dMinor(end, :, 2), 'r');
handles.PlotdMinorl = plot(xAB(1,:), dMinor(end, :, 1), 'b');
xlim([0 xAB(1,end)]);
ylim([0 max(dMinor(:))]);
handles.PlotDisbondPos1 = plot([b(end) b(end)], ylim, '-mo');
hold off
title('Minor Index Increment (\DeltaMinor / \Deltab)')
xlabel('x_{AB}')
ylabel('\DeltaMinor / \Deltab [-]')
legend('Top Ply', 'Bottom ply', 'Disbond [mm]')
grid on

% Update handles structure
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function axes6_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function axes3_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function axes1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
