function varargout = menu_process(varargin)
% MENU_PROCESS M-file for menu_process.fig
%      MENU_PROCESS, by itself, creates a new MENU_PROCESS or raises the existing
%      singleton*.
%
%      H = MENU_PROCESS returns the handle to a new MENU_PROCESS or the handle to
%      the existing singleton*.
%
%      MENU_PROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MENU_PROCESS.M with the given input arguments.
%
%      MENU_PROCESS('Property','Value',...) creates a new MENU_PROCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before menu_process_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to menu_process_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help menu_process

% Last Modified by GUIDE v2.5 28-Feb-2007 09:26:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @menu_process_OpeningFcn, ...
    'gui_OutputFcn',  @menu_process_OutputFcn, ...
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
set(handles.browserf, 'UserData', handles.editloadfrom);
set(handles.browserp, 'UserData', handles.editsaveto);
step = getappdata(0,'step');
frequency = getappdata(0,'frequency');
begin = getappdata(0,'begin');
ndp = getappdata(0,'ndp');
[allow_damp,allow_freq,lambda,maxiter,phasedistort,filterlength] = computeAQSESparam(step,begin,ndp,frequency);
set(handles.editfilterlength,'String',filterlength);
set(handles.editallowdamp,'String',allow_damp);
set(handles.editallowfreq,'String',allow_freq);
set(handles.editlambda,'String',lambda);
set(handles.editmaxiter,'String',maxiter);
set(handles.editECC,'String',phasedistort);
% --- Executes just before menu_process is made visible.
function menu_process_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to menu_process (see VARARGIN)

% Choose default command line output for menu_process
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes menu_process wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = menu_process_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in processbutton.
function processbutton_Callback(hObject, eventdata, handles)
% hObject    handle to processbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ndp = getappdata(0,'ndp');
signal = getappdata(0,'signal');
step = getappdata(0,'step');
begin = getappdata(0,'begin');
frequency = getappdata(0,'frequency');
if strcmp(frequency,'')
    frequency = 63130;
    warning(['The frequency has been defined as ',num2str(frequency),' kHz.'])
end

typemethode = get(handles.quantmethod,'Value'); %typemethode = 1 if AQSES
%         = 2 if HLSVD-PRO

filtertype = get(handles.editfiltertype,'Value')-1; %  0 = no filter;filtertype = 1 if MP-FIR
%        = 2 if HLSVD-PRO, =3 if MP-FIR0 by (JB Poullet)


% if  filtertype==1
%     filtertype = '';
% elseif filtertype==2
%     filtertype = 'fir';
% elseif filtertype==3
%     filtertype = 'hlsvd';
% end

allow_damp = str2num(get(handles.editallowdamp,'String'));
allow_freq = str2num(get(handles.editallowfreq,'String'))*2*pi;
lambda = str2num(get(handles.editlambda,'String'));
maxiter = str2num(get(handles.editmaxiter,'String'));
phasedistort = str2num(get(handles.editECC,'String'));


ppmrange = str2num(get(handles.editrange,'String'));
boundL = ppmrange(1);
boundH = ppmrange(2);
fsampl = 1/step;
%xmin_mouse = ((ppmrange(1)-4.7)/(fsampl/(frequency/1000)))/1000;%in kHz
%xmax_mouse = ((ppmrange(2)-4.7)/(fsampl/(frequency/1000)))/1000;%in kHz
ripple = str2num(get(handles.editripple,'String'));
filterlength = str2num(get(handles.editfilterlength,'String'));
equalph = get(handles.checkequalph,'Value');
baseline_boolean = get(handles.checkbaseline,'Value');

if get(handles.checkplotiter,'Value')
    plotting = 1;
else
    plotting = 0;
end


truncation = get(handles.checktruncation,'Value');
if truncation
    trunc = str2num(get(handles.edittrunc,'String'));
    st = get(handles.beginend,'String');
    truncationtype = st(get(handles.beginend,'Value'));
else
    trunc = 0;
    truncationtype =  'Begin';
end

if truncation
    signal = signal(1:length(signal)-trunc);
    ndp = ndp-trunc;
end

t   = [0:step:(ndp-1)*step] + begin;%in ms

% construct database
% ndpt = ndp;
% begint = begin;
% stept = step;


ff  = get(handles.editloadfrom,'String');
dbase = get(handles.editloadfrom,'String');
if strcmp(ff,'')
    errordlg('You must enter a database file');
    return
else

    abs_scl_factor = 1;
    lineshape = 1;
    prior_equal = 0;
    equal_to = 0;

    %NEW PART
    for i =  1:size(signal,1)
        [scores(i,:),misc(i)] = AQSES(signal(i,:),step,frequency,ndp,begin,dbase,abs_scl_factor,...
            lineshape, phasedistort,filtertype,boundL,boundH,ripple,filterlength,...
            prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
            maxiter,plotting);
    end

    procres.scores = scores;
    procres.misc = misc;
    classopt = [];
    savef = get(handles.editsaveto,'String');
    if strcmp(savef,'')
        errordlg('You must first specify a file in which you want to store the results.')
    else
        save(savef,'procres','classopt','signal','ndp','begin','step','frequency');
    end
end


% --- Executes on selection change in quantmethod.
function quantmethod_Callback(hObject, eventdata, handles)
% hObject    handle to quantmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns quantmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from quantmethod


% --- Executes during object creation, after setting all properties.
function quantmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quantmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function editallowfreq_Callback(hObject, eventdata, handles)
% hObject    handle to editallowfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editallowfreq as text
%        str2double(get(hObject,'String')) returns contents of editallowfreq as a double


% --- Executes during object creation, after setting all properties.
function editallowfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editallowfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editlambda_Callback(hObject, eventdata, handles)
% hObject    handle to editlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlambda as text
%        str2double(get(hObject,'String')) returns contents of editlambda as a double


% --- Executes during object creation, after setting all properties.
function editlambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function editmaxiter_Callback(hObject, eventdata, handles)
% hObject    handle to editmaxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmaxiter as text
%        str2double(get(hObject,'String')) returns contents of editmaxiter as a double


% --- Executes during object creation, after setting all properties.
function editmaxiter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmaxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editECC_Callback(hObject, eventdata, handles)
% hObject    handle to editECC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editECC as text
%        str2double(get(hObject,'String')) returns contents of editECC as a double


% --- Executes during object creation, after setting all properties.
function editECC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editECC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function editrange_Callback(hObject, eventdata, handles)
% hObject    handle to editrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editrange as text
%        str2double(get(hObject,'String')) returns contents of editrange as a double


% --- Executes during object creation, after setting all properties.
function editrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editripple_Callback(hObject, eventdata, handles)
% hObject    handle to editripple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editripple as text
%        str2double(get(hObject,'String')) returns contents of editripple as a double


% --- Executes during object creation, after setting all properties.
function editripple_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editripple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in editfiltertype.
function editfiltertype_Callback(hObject, eventdata, handles)
% hObject    handle to editfiltertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns editfiltertype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from editfiltertype


% --- Executes during object creation, after setting all properties.
function editfiltertype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editfiltertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editfilterlength_Callback(hObject, eventdata, handles)
% hObject    handle to editfilterlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editfilterlength as text
%        str2double(get(hObject,'String')) returns contents of editfilterlength as a double


% --- Executes during object creation, after setting all properties.
function editfilterlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editfilterlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on selection change in beginend.
function beginend_Callback(hObject, eventdata, handles)
% hObject    handle to beginend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns beginend contents as cell array
%        contents{get(hObject,'Value')} returns selected item from beginend


% --- Executes during object creation, after setting all properties.
function beginend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beginend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edittrunc_Callback(hObject, eventdata, handles)
% hObject    handle to edittrunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittrunc as text
%        str2double(get(hObject,'String')) returns contents of edittrunc as a double


% --- Executes during object creation, after setting all properties.
function edittrunc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittrunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in checktruncation.
function checktruncation_Callback(hObject, eventdata, handles)
% hObject    handle to checktruncation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checktruncation





function editloadfrom_Callback(hObject, eventdata, handles)
% hObject    handle to editloadfrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editloadfrom as text
%        str2double(get(hObject,'String')) returns contents of editloadfrom as a double
loadmetabolites(hObject, eventdata, handles);


function loadmetabolites(hObject, eventdata, handles)
set(handles.editmetabolites,'String','');
ff  = get(handles.editloadfrom,'String');
a = load(ff);
aN = fieldnames(a);
c=0;

if isfield(a,'ndp')
    for i=1:length(aN)
        if ~strcmp(char(aN(i)),'step')
            if size(eval(['a.',char(aN(i))]),2)==a.ndp
                cem = get(handles.editmetabolites,'String');
                set(handles.editmetabolites,'String',[cem, ' ',char(aN(i))]);
                %set(handles.editmetplot,'String',[cem, ' ',char(aN(i))]);
                c=c+1;
            end
        end
    end
else
    errordlg('The database should contain at least one metabolite profile and the variable ndp');
end
set(handles.editnbrmet,'String',num2str(c));




% --- Executes during object creation, after setting all properties.
function editloadfrom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editloadfrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in browserf.
function browserf_Callback(hObject, eventdata, handles)
% hObject    handle to browserf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=get(gcbo,'UserData');
[fn p]=uigetfile(get(f, 'String'),'Load file from');
if isequal(fn,0) | isequal(p,0)
    return
end
cd(pwd);
if (~isempty(fn)&fn~=0),
    set(f,'String',fullfile(p, fn));
end
loadmetabolites(hObject, eventdata, handles);


% --- Executes on button press in checkbaseline.
function checkbaseline_Callback(hObject, eventdata, handles)
% hObject    handle to checkbaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbaseline





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




function editmetabolites_Callback(hObject, eventdata, handles)
% hObject    handle to editmetabolites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmetabolites as text
%        str2double(get(hObject,'String')) returns contents of editmetabolites as a double
[A, count, errmsg, nextindex] = sscanf(get(handles.editmetabolites,'String'),'%s');
set(handles.editnbrmet,'String',count)

% --- Executes during object creation, after setting all properties.
function editmetabolites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmetabolites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function editnbrmet_Callback(hObject, eventdata, handles)
% hObject    handle to editnbrmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editnbrmet as text
%        str2double(get(hObject,'String')) returns contents of editnbrmet as a double


% --- Executes during object creation, after setting all properties.
function editnbrmet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editnbrmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in checkplotiter.
function checkplotiter_Callback(hObject, eventdata, handles)
% hObject    handle to checkplotiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkplotiter




% --- Executes on button press in checkmetplot.
function checkmetplot_Callback(hObject, eventdata, handles)
% hObject    handle to checkmetplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkmetplot


function editmetplot_Callback(hObject, eventdata, handles)
% hObject    handle to editmetplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmetplot as text
%        str2double(get(hObject,'String')) returns contents of editmetplot as a double


% --- Executes during object creation, after setting all properties.
function editmetplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmetplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in checkores.
function checkores_Callback(hObject, eventdata, handles)
% hObject    handle to checkores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkores





function editallowdamp_Callback(hObject, eventdata, handles)
% hObject    handle to editallowdamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editallowdamp as text
%        str2double(get(hObject,'String')) returns contents of editallowdamp as a double


% --- Executes during object creation, after setting all properties.
function editallowdamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editallowdamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in checkequalph.
function checkequalph_Callback(hObject, eventdata, handles)
% hObject    handle to checkequalph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkequalph


