function varargout = GUI_FinalPlots(varargin)
% GUI_FINALPLOTS MATLAB code for GUI_FinalPlots.fig
%      GUI_FINALPLOTS, by itself, creates a new GUI_FINALPLOTS or raises the existing
%      singleton*.
%
%      H = GUI_FINALPLOTS returns the handle to a new GUI_FINALPLOTS or the handle to
%      the existing singleton*.
%
%      GUI_FINALPLOTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FINALPLOTS.M with the given input arguments.
%
%      GUI_FINALPLOTS('Property','Value',...) creates a new GUI_FINALPLOTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FinalPlots_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FinalPlots_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FinalPlots

% Last Modified by GUIDE v2.5 23-May-2018 16:48:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FinalPlots_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FinalPlots_OutputFcn, ...
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


% --- Executes just before GUI_FinalPlots is made visible.
function GUI_FinalPlots_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FinalPlots (see VARARGIN)

% Choose default command line output for GUI_FinalPlots
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_FinalPlots wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FinalPlots_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function Slider_Element_Callback(hObject, eventdata, handles)
% hObject    handle to Slider_Element (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


Shear_a = getappdata(0,'Shear_a')*1e-6;
Peel_a  = getappdata(0,'Peel_a')*1e-6;
N1      = getappdata(0,'N1')*1e-3;
N2      = getappdata(0,'N2')*1e-3;
M1      = getappdata(0,'M1')*1e-3;
M2      = getappdata(0,'M2')*1e-3;
Q1      = getappdata(0,'Q1')*1e-3;
Q2      = getappdata(0,'Q2')*1e-3;
x       = getappdata(0, 'x');

value = get(hObject, 'value');
value = round(value);

% Update the plots
set(handles.PlotShear1,'XData',x(value, value:end),'YData',Shear_a(value, value:end,1))
set(handles.PlotShear2,'XData',x(value, value:end),'YData',Shear_a(value, value:end,2))

set(handles.PlotPeel1,'XData',x(value, value:end),'YData',Peel_a(value, value:end,1))
set(handles.PlotPeel2,'XData',x(value, value:end),'YData',Peel_a(value, value:end,2))

set(handles.PlotM1,'XData',x(value, value:end),'YData',M1(value, value:end,2))
set(handles.PlotM2,'XData',x(value, value:end),'YData',M2(value, value:end,2))

set(handles.PlotN1,'XData',x(value, value:end),'YData',N1(value, value:end,2))
set(handles.PlotN2,'XData',x(value, value:end),'YData',N2(value, value:end,2))

set(handles.PlotQ1,'XData',x(value, value:end),'YData',Q1(value, value:end,2))
set(handles.PlotQ2,'XData',x(value, value:end),'YData',Q2(value, value:end,2))

% --- Executes during object creation, after setting all properties.
function Slider_Element_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slider_Element (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

if isappdata(0, 'x')
    x = getappdata(0, 'x');
else
    x = evalin('base','x00');
    setappdata(0, 'x', x);
end
    
set(hObject, 'Min', 1);
set(hObject, 'Max', size(x,1));
set(hObject, 'Value', 1);
set(hObject, 'SliderStep', [10/(size(x,1)-1) , 100/(size(x,1)-1)]);

% --- Executes during object creation, after setting all properties.
function axes_AdhesiveShear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_AdhesiveShear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_AdhesiveShear
% if isappdata(0, 'Shear_a')
%     Shear_a = getappdata(0, 'Shear_a');
% else
%     Shear_a = evalin('base','Shear_a');
%     setappdata(0, 'Shear_a', Shear_a);
% end
% 
% if isappdata(0, 'x')
%     x = getappdata(0, 'x');
% else
%     x = evalin('base','x00');
%     setappdata(0, 'x', x);
% end

Shear_a = getappdata(0, 'Shear_a')*1e-6;
x = getappdata(0, 'x');

q = 1;

hold on
handles.PlotShear1 = plot(x(q,:),Shear_a(q,:,1),'b');
handles.PlotShear2 = plot(x(q,:),Shear_a(q,:,2),'r');
hold off
title('Adhesive Shear Stress (\tau_{xy})')
ylabel('\tau_{xy} [MPa]')
xlabel('Distance from the right overlap edge [mm]')
legend('\tau_{xy} at S_{min}', '\tau_{xy} at S_{max}')
xlim([x(1,1) x(1,end)]);
ylim([min(Shear_a(:)) max(Shear_a(:))]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_AdhesivePeel_CreateFcn(hObject, ~, handles)
% hObject    handle to axes_AdhesivePeel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_AdhesivePeel
% if isappdata(0, 'Peel_a')
%     Peel_a = getappdata(0, 'Peel_a');
% else
%     Peel_a = evalin('base','Peel_a');
%     setappdata(0, 'Peel_a', Peel_a);
% end
% 
% if isappdata(0, 'x')
%     x = getappdata(0, 'x');
% else
%     x = evalin('base','x00');
%     setappdata(0, 'x', x);
% end

Peel_a = getappdata(0, 'Peel_a')*1e-6;
x = getappdata(0, 'x');

q = 1;

hold on
handles.PlotPeel1 = plot(x(q,:),Peel_a(q,:,1),'b');
handles.PlotPeel2 = plot(x(q,:),Peel_a(q,:,2),'r');
hold off
title('Adhesive Peel Stress (\sigma_{y})')
ylabel('\sigma_{y} [MPa]')
xlabel('Distance from the right overlap edge [mm]')
legend('\sigma_{y} at S_{min}', '\sigma_{y} at S_{max}')
xlim([x(1,1) x(1,end)]);
ylim([min(Peel_a(:)) max(Peel_a(:))]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_M_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_M
% if isappdata(0, 'M1')
%     M1 = getappdata(0, 'M1');
% else
%     M1 = evalin('base','Loads_ad1.M');
%     setappdata(0, 'M1', M1);
% end
% 
% if isappdata(0, 'M2')
%     M2 = getappdata(0, 'M2');
% else
%     M2 = evalin('base','Loads_ad2.M');
%     setappdata(0, 'M2', M2);
% end
% 
% if isappdata(0, 'x')
%     x = getappdata(0, 'x');
% else
%     x = evalin('base','x00');
%     setappdata(0, 'x', x);
% end

M1 = getappdata(0, 'M1')*1e-3;
M2 = getappdata(0, 'M2')*1e-3;
x = getappdata(0, 'x');

q = 1;

hold on
handles.PlotM1 = plot(x(q,:),M1(q,:,2),'b');
handles.PlotM2 = plot(x(q,:),M2(q,:,2),'r');
hold off
title('Adherent Bending Moment (M) at S_{max}')
ylabel('M [kNm]')
xlabel('Distance from the right overlap edge [mm]')
legend('Upper adherent', 'Lower adherent')
xlim([x(1,1) x(1,end)]);
ylim([min(M2(:)) max(M1(:))]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_N
% if isappdata(0, 'N1')
%     N1 = getappdata(0, 'Q1');
% else
%     N1 = evalin('base','Loads_ad1.N');
%     setappdata(0, 'N1', N1);
% end
% 
% if isappdata(0, 'N2')
%     N2 = getappdata(0, 'N2');
% else
%     N2 = evalin('base','Loads_ad2.N');
%     setappdata(0, 'N2', N2);
% end
% 
% if isappdata(0, 'x')
%     x = getappdata(0, 'x');
% else
%     x = evalin('base','x00');
%     setappdata(0, 'x', x);
% end

N1 = getappdata(0, 'N1')*1e-3;
N2 = getappdata(0, 'N2')*1e-3;
x = getappdata(0, 'x');

q = 1;

hold on
handles.PlotN1 = plot(x(q,:),N1(q,:,2),'b');
handles.PlotN2 = plot(x(q,:),N2(q,:,2),'r');
hold off
title('Adherent Axial Force (N) at S_{max}')
ylabel('N [kN/m]')
xlabel('Distance from the right overlap edge [mm]')
legend('Upper adherent', 'Lower adherent')
xlim([x(1,1) x(1,end)]);
ylim([min(N2(:)) max(N1(:))]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_Q
% if isappdata(0, 'Q1')
%     Q1 = getappdata(0, 'Q1');
% else
%     Q1 = evalin('base','Loads_ad1.Q');
%     setappdata(0, 'Q1', Q1);
% end
% 
% if isappdata(0, 'Q2')
%     Q2 = getappdata(0, 'Q2');
% else
%     Q2 = evalin('base','Loads_ad2.Q');
%     setappdata(0, 'Q2', Q2);
% end
% 
% if isappdata(0, 'x')
%     x = getappdata(0, 'x');
% else
%     x = evalin('base','x00');
%     setappdata(0, 'x', x);
% end

Q1 = getappdata(0, 'Q1')*1e-3;
Q2 = getappdata(0, 'Q2')*1e-3;
x = getappdata(0, 'x');

q = 1;

hold on
handles.PlotQ1 = plot(x(q,:),Q1(q,:,2),'b');
handles.PlotQ2 = plot(x(q,:),Q2(q,:,2),'r');
hold off
title('Adherent Shear Force (Q) at S_{max}')
ylabel('Q [kN/m]')
xlabel('Distance from the right overlap edge [mm]')
legend('Upper adherent', 'Lower adherent')
xlim([x(1,1) x(1,end)]);
ylim([min(Q1(:)) max(Q2(:))]);
grid on

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes_CLS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_CLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_CLS



% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
