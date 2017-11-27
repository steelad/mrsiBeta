function varargout = menu_peakintegration(varargin)
% MENU_PEAKINTEGRATION M-file for menu_peakintegration.fig
%      MENU_PEAKINTEGRATION, by itself, creates a new MENU_PEAKINTEGRATION or raises the existing
%      singleton*.
%
%      H = MENU_PEAKINTEGRATION returns the handle to a new MENU_PEAKINTEGRATION or the handle to
%      the existing singleton*.
%
%      MENU_PEAKINTEGRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MENU_PEAKINTEGRATION.M with the given input arguments.
%
%      MENU_PEAKINTEGRATION('Property','Value',...) creates a new MENU_PEAKINTEGRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before menu_peakintegration_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to menu_peakintegration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help menu_peakintegration

% Last Modified by GUIDE v2.5 02-Mar-2007 11:11:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @menu_peakintegration_OpeningFcn, ...
    'gui_OutputFcn',  @menu_peakintegration_OutputFcn, ...
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


% --- Executes just before menu_peakintegration is made visible.
function menu_peakintegration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to menu_peakintegration (see VARARGIN)

% Choose default command line output for menu_peakintegration
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes menu_peakintegration wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(hObject, handles);
end


function initialize_gui(fig_handle, handles)
set(handles.browserp, 'UserData', handles.editsaveto);
[compNames,freqmatrix] = readintervals(1);
dispmatrix(freqmatrix,0)
cem = [];
for i=1:length(compNames)
    cem = [cem, ' ',char(compNames(i))];
end
set(handles.editcompNames,'String',cem)
%set(handles.editintervals,'Value',freqmatrix)

% --- Outputs from this function are returned to the command line.
function varargout = menu_peakintegration_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in integ_range.
function integ_range_Callback(hObject, eventdata, handles)
% hObject    handle to integ_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns integ_range contents as cell array
%        contents{get(hObject,'Value')} returns selected item from integ_range
val = get(handles.integ_range,'Value');
str = get(handles.integ_range,'String');

if strcmp(char(str(val)),'Other')
    set(handles.editcompNames,'String','{''met1'' ''met2''}')
    set(handles.editintervals,'String','[0 0 0 0]')
    dispmatrix([0 0; 0 0],0);
else
    [compNames,freqmatrix] = readintervals(val);
    dispmatrix(freqmatrix,0)
    cem = [];
    for i=1:length(compNames)
        cem = [cem, ' ',char(compNames(i))];
    end
    set(handles.editcompNames,'String',cem)
    cem = [];
    for i=1:size(freqmatrix,1)
        cem = [cem, '   ',num2str(freqmatrix(i,1)),' ',num2str(freqmatrix(i,2))];
    end
    set(handles.editintervals,'String',cem)
end


% --- Executes during object creation, after setting all properties.
function integ_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to integ_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editcompNames_Callback(hObject, eventdata, handles)
% hObject    handle to editcompNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editcompNames as text
%        str2double(get(hObject,'String')) returns contents of editcompNames as a double


% --- Executes during object creation, after setting all properties.
function editcompNames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editcompNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in peakintegration.
function peakintegration_Callback(hObject, eventdata, handles)
% hObject    handle to peakintegration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
step = getappdata(0,'step');
signal = getappdata(0,'signal');
frequency = getappdata(0,'frequency');
begin = getappdata(0,'begin');
ndp = size(signal,2);

val = get(handles.integ_range,'Value')
str = get(handles.integ_range,'String')
freqmatrix = str2num(get(handles.editintervals,'String'))
compNames = cellstr(get(handles.editcompNames,'String'))

if strcmp(char(str(val)),'Other')
    [scores,misc] = peakintegtempl(signal,step,frequency,ndp,begin,val,freqmatrix,compNames);
else
    [scores,misc] = peakintegtempl(signal,step,frequency,ndp,begin,val);
end
procres.scores = scores;
procres.misc = misc;
savef = get(handles.editsaveto,'String');
save(savef,'procres');



function editintervals_Callback(hObject, eventdata, handles)
% hObject    handle to editintervals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editintervals as text
%        str2double(get(hObject,'String')) returns contents of editintervals as a double
freqmatrix = str2num(get(handles.editintervals,'string'));
for i = 1:length(freqmatrix)/2
    freqmatrixt(i,1:2) = freqmatrix((i-1)*2+1:i*2);
end
freqmatrix = freqmatrixt;
dispmatrix(freqmatrix,0);

% --- Executes during object creation, after setting all properties.
function editintervals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editintervals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function editsaveto_Callback(hObject, eventdata, handles)
% hObject    handle to editsaveto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editsaveto as text
%        str2double(get(hObject,'String')) returns contents of editsaveto as a double


% --- Executes during object creation, after setting all properties.
function editsaveto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editsaveto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in browserp.
function browserp_Callback(hObject, eventdata, handles)
% hObject    handle to browserp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=get(gcbo,'UserData');
[fn p]=uiputfile(get(f, 'String'), 'Save to');
if isequal(fn,0) | isequal(p,0)
    return
end
cd(pwd);
if (~isempty(fn)&fn~=0),
    set(f,'String',fullfile(p, fn));
end

