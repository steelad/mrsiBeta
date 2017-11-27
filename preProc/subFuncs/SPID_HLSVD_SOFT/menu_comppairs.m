function varargout = menu_compmethods(varargin)
% MENU_COMPMETHODS M-file for menu_compmethods.fig
%      MENU_COMPMETHODS, by itself, creates a new MENU_COMPMETHODS or raises the existing
%      singleton*.
%
%      H = MENU_COMPMETHODS returns the handle to a new MENU_COMPMETHODS or the handle to
%      the existing singleton*.
%
%      MENU_COMPMETHODS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MENU_COMPMETHODS.M with the given input arguments.
%
%      MENU_COMPMETHODS('Property','Value',...) creates a new MENU_COMPMETHODS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before menu_compmethods_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to menu_compmethods_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help menu_compmethods

% Last Modified by GUIDE v2.5 30-Oct-2007 17:01:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @menu_comppairs_OpeningFcn, ...
                   'gui_OutputFcn',  @menu_comppairs_OutputFcn, ...
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
function initialize_gui(fig_handle, handles)
set(handles.browserp, 'UserData', handles.editsaveto);

% --- Executes just before menu_compmethods is made visible.
function menu_compmethods_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to menu_compmethods (see VARARGIN)

% Choose default command line output for menu_compmethods
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes menu_compmethods wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = menu_compmethods_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in comppair.
function comppair_Callback(hObject, eventdata, handles)
% hObject    handle to comppair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=get(gcbo,'UserData');
[fn p]=uigetfile(get(f, 'String'),'Load file from');
if isequal(fn,0) | isequal(p,0)
    return
end
ttestcomp(fullfile(p,char(fn)));

% --- Executes on selection change in compmethods.
function compmethods_Callback(hObject, eventdata, handles)
% hObject    handle to compmethods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns compmethods contents as cell array
%        contents{get(hObject,'Value')} returns selected item from compmethods


% --- Executes during object creation, after setting all properties.
function compmethods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compmethods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


