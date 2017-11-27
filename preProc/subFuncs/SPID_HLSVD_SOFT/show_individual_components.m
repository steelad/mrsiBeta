function varargout = show_individual_components(varargin)
% SHOW_INDIVIDUAL_COMPONENTS M-file for show_individual_components.fig
%      SHOW_INDIVIDUAL_COMPONENTS, by itself, creates a new SHOW_INDIVIDUAL_COMPONENTS or raises the existing
%      singleton*.
%
%      H = SHOW_INDIVIDUAL_COMPONENTS returns the handle to a new SHOW_INDIVIDUAL_COMPONENTS or the handle to
%      the existing singleton*.
%
%      SHOW_INDIVIDUAL_COMPONENTS('Property','Value',...) creates a new SHOW_INDIVIDUAL_COMPONENTS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to show_individual_components_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SHOW_INDIVIDUAL_COMPONENTS('CALLBACK') and SHOW_INDIVIDUAL_COMPONENTS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SHOW_INDIVIDUAL_COMPONENTS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_individual_components

% Last Modified by GUIDE v2.5 22-Jul-2004 18:03:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_individual_components_OpeningFcn, ...
                   'gui_OutputFcn',  @show_individual_components_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before show_individual_components is made visible.
function show_individual_components_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for show_individual_components
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show_individual_components wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show_individual_components_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in names.
function names_Callback(hObject, eventdata, handles, ...
    ampl_ratios,ampl_absolute,CR_error,individual_components,...
    phas,step, t0, fir_h, xmin_mouse, xmax_mouse)
% hObject    handle to names (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns names contents as cell array
%        contents{get(hObject,'Value')} returns selected item from names

which_metabolite = get(hObject,'Value');
set(handles.ratios,'String',ampl_ratios(which_metabolite))
set(handles.concentrations,'String',ampl_absolute(which_metabolite))
set(handles.errors,'String',CR_error(which_metabolite))
set(handles.disp_components,'Visible','on')
ndft = size(individual_components,2);
axes(handles.disp_components)
disp_sigreal(individual_components(which_metabolite,:).',phas,step,t0,ndft,' ',0,3,fir_h); 
v = axis;
xmin_mouse=4.65+1000*xmin_mouse/63.86;
xmax_mouse=4.65+1000*xmax_mouse/63.86;
axis([xmin_mouse xmax_mouse v(3) v(4)]);  
grid, hold off 

%%%% must be set in the ButtonDwn callback:
%%%% show_individual_components('names_Callback',gcbo,[],guidata(gcbo),ampl_ratios,ampl_absolute,CR_error,individual_components,step, t0, fir_h, xmin_mouse, xmax_mouse)


% --- Executes during object creation, after setting all properties.
function ratios_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ratios (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'String','')

% --- Executes during object creation, after setting all properties.
function concentrations_CreateFcn(hObject, eventdata, handles)

set(hObject,'String','')

% --- Executes during object creation, after setting all properties.
function errors_CreateFcn(hObject, eventdata, handles)

set(hObject,'String','')

