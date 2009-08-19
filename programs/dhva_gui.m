function varargout = dHvA_GUI(varargin)
% DHVA_GUI M-file for dHvA_GUI.fig
%      DHVA_GUI, by itself, creates a new DHVA_GUI or raises the existing
%      singleton*.
%
%      H = DHVA_GUI returns the handle to a new DHVA_GUI or the handle to
%      the existing singleton*.
%
%      DHVA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DHVA_GUI.M with the given input arguments.
%
%      DHVA_GUI('Property','Value',...) creates a new DHVA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dHvA_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dHvA_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help dHvA_GUI

% Last Modified by GUIDE v2.5 14-Nov-2005 17:36:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dHvA_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @dHvA_GUI_OutputFcn, ...
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

% --- Executes just before dHvA_GUI is made visible.
function dHvA_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dHvA_GUI (see VARARGIN)

% Choose default command line output for dHvA_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dHvA_GUI wait for user response (see UIRESUME)
% uiwait(handles.MainWindow);


% --- Outputs from this function are returned to the command line.
function varargout = dHvA_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function CloseMainWindow_Callback(hObject, eventdata, handles)
global actions;
actions='Close';

function Axes_Rotplot_ButtonDown_Callback(hObject, eventdata, handles)
global actions;
actions='Clicked rotplot';

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MI_OpenOrbitsFile_Callback(hObject, eventdata, handles)
% hObject    handle to MI_OpenOrbitsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global actions;

actions='Open orbits file';

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.MainWindow)

function CB_BackProject_Callback(hObject, eventdata, handles)
global actions;

actions='Changed backproj';

function ET_MagField_Callback(hObject, eventdata, handles)
global actions;

%actions='Changed backproj';

% --------------------------------------------------------------------
function MI_Close_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global actions;
actions='Close';

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on selection change in LB_bandsdatafiles.
function LB_bandsdatafiles_Callback(hObject, eventdata, handles)
% hObject    handle to LB_bandsdatafiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LB_bandsdatafiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LB_bandsdatafiles
global actions;
actions='Selected bandsdatafile';

% --- Executes during object creation, after setting all properties.
function LB_bandsdatafiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_bandsdatafiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LB_bands.
function LB_bands_Callback(hObject, eventdata, handles)
% hObject    handle to LB_bands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LB_bands contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LB_bands
global actions;
actions='Selected band';

% --- Executes during object creation, after setting all properties.
function LB_bands_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_bands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BG_Undersampling_SelectionChangeFcn(hObject,eventdata,handles)
% hObject    handle to uipanel1 (see GCBO) WRONG!
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global actions;
allbuttons=get(get(hObject,'Parent'),'Children');
set(allbuttons,'Value',0);
set(hObject,'Value',1);
actions='Changed undersampling';

function BG_PlotRegion_SelectionChangeFcn(hObject,eventdata,handles)
% hObject    handle to uipanel1 (see GCBO) WRONG!
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global actions;
allbuttons=get(get(hObject,'Parent'),'Children');
set(allbuttons,'Value',0);
set(hObject,'Value',1);
actions='Changed clipping';

% --- Executes on selection change in LB_ExtremalOrbits.
function LB_ExtremalOrbits_Callback(hObject, eventdata, handles)
% hObject    handle to LB_ExtremalOrbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LB_ExtremalOrbits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LB_ExtremalOrbits
global actions;
actions='Selected extremal orbit';

% --- Executes during object creation, after setting all properties.
function LB_ExtremalOrbits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_ExtremalOrbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in LB_orbitsdatafiles.
function LB_orbitsdatafiles_Callback(hObject, eventdata, handles)
% hObject    handle to LB_orbitsdatafiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LB_orbitsdatafiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LB_orbitsdatafiles
global actions;
actions='Selected orbitsdatafile';

% --- Executes during object creation, after setting all properties.
function LB_orbitsdatafiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LB_orbitsdatafiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function MI_Expand_klist_Callback(hObject, eventdata, handles)
% hObject    handle to MI_Expand_klist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MI_Convert_Wien2k_Callback(hObject, eventdata, handles)
% hObject    handle to MI_Convert_Wien2k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MI_CalcRotPlot_Callback(hObject, eventdata, handles)
% hObject    handle to MI_CalcRotPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes when MainWindow is resized.
function MainWindow_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to MainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ET_dMdB_Callback(hObject, eventdata, handles)

global actions;

actions='Changed backproj';


% --- Executes during object creation, after setting all properties.
function ET_dMdB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_dMdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ET_MagField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_MagField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CB_ShowBZ.
function CB_ShowBZ_Callback(hObject, eventdata, handles)
% hObject    handle to CB_ShowBZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_ShowBZ
global actions;
actions='Clicked ShowBZ';



% --- Executes on button press in CB_ShowContourPlot.
function CB_ShowContourPlot_Callback(hObject, eventdata, handles)
% hObject    handle to CB_ShowContourPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_ShowContourPlot
global actions;
actions='Clicked ShowContourplot';




function ET_ContourPlotBands_Callback(hObject, eventdata, handles)
% hObject    handle to ET_ContourPlotBands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_ContourPlotBands as text
%        str2double(get(hObject,'String')) returns contents of ET_ContourPlotBands as a double


% --- Executes during object creation, after setting all properties.
function ET_ContourPlotBands_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_ContourPlotBands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


