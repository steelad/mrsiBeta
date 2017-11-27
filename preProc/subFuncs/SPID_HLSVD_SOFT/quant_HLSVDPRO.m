function varargout = quant_HLSVDPRO(varargin)
% QUANT_HLSVDPRO M-file for quant_HLSVDPRO.fig
%      QUANT_HLSVDPRO, by itself, creates a new QUANT_HLSVDPRO or raises the existing
%      singleton*.
%
%      H = QUANT_HLSVDPRO returns the handle to a new QUANT_HLSVDPRO or the handle to
%      the existing singleton*.
%
%      QUANT_HLSVDPRO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUANT_HLSVDPRO.M with the given input arguments.
%
%      QUANT_HLSVDPRO('Property','Value',...) creates a new QUANT_HLSVDPRO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quant_HLSVDPRO_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quant_HLSVDPRO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help quant_HLSVDPRO

% Last Modified by GUIDE v2.5 23-Nov-2006 17:19:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @quant_HLSVDPRO_OpeningFcn, ...
    'gui_OutputFcn',  @quant_HLSVDPRO_OutputFcn, ...
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


% --- Executes just before quant_HLSVDPRO is made visible.
function quant_HLSVDPRO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quant_HLSVDPRO (see VARARGIN)

% Choose default command line output for quant_HLSVDPRO
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes quant_HLSVDPRO wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(hObject, handles);
end

function initialize_gui(fig_handle, handles)
set(handles.browserp, 'UserData', handles.editsaveto);



% --- Outputs from this function are returned to the command line.
function varargout = quant_HLSVDPRO_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in quant.
function quant_Callback(hObject, eventdata, handles)
% hObject    handle to quant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ndp = getappdata(0,'ndp');
signal = getappdata(0,'signal');
step = getappdata(0,'step');
begin = getappdata(0,'begin');
frequency = getappdata(0,'frequency');
signalold = signal;
nos = size(signal,1);
if strcmp(frequency,'')
    frequency = 63130;
    warning(['The frequency has been defined as ',num2str(frequency),' kHz.'])
end

ppmrange = str2num(get(handles.editrange,'String'));
fsampl = 1/step;
boundL = ppmrange(1);
boundH = ppmrange(2);
% xmin = (ppmrange(1)-4.7)*frequency/10^6;%in kHz
% xmax = (ppmrange(2)-4.7)*frequency/10^6;%in kHz
Kest = str2num(get(handles.editorder,'String'));

amplAll = []; freqAll = []; phasAll = []; dampAll=[];
filoriginalAll = []; filbaselineAll = []; filreconAll=[];

if get(handles.method,'value')==1 %HLSVDPRO
    
    for i=1:nos
        %t = [0:step:(ndp-1)*step];
        [scores(i,:),misc(i)] = HLSVDPROquant(signal(i,:),step,frequency,ndp,begin,boundL,boundH,Kest);

        %         freqrange1 = [-(fsampl*1000/2-1),xmin*1000]; %in Hz
        %         freqrange2 = [xmax*1000,(fsampl*1000/2-1)]; %in Hz
        %         freqrange = [freqrange1;freqrange2];
        %         %     freqrange1 = [-499,-280]; %in Hz
        %         %     freqrange2 = [-100,499]; %in Hz
        %         MM = ndp/2;
        %         fsampl = 1/step*1000; %in Hz
        %         [freq,damp,ampl,phas]=HLSVDPRO(signal(i,:),signal(i,:),Kest,fsampl,t/1000,MM);
        %           [filrecon(i,:)]=reconsd(t/1000,misc(i).freq,misc(i).damp,misc(i).ampl,misc(i).phas); %reconstructed signals
        %           misc(i).filoriginal = signal(i,:); % these signals are the original ones
        %         freqAll(i,:) = freq/1000; %to get the frequencies in kHz
        %         amplAll(i,:) = ampl;
        %         dampAll(i,:) = damp/1000;
        %         phasAll(i,:) = phas;
        %         filreconAll(i,:) = filrecon;
        %         for k = 1:Kest
        %             individual_components(k,:) = reconsd(t/1000,misc(i).freq(k),misc(i).damp(k),misc(i).ampl(k),misc(i).phas(k));
        %         end
        %         individual_componentsAll(i).indcomp = individual_components;
    end %i=1:nos
    compNames(1,:) =['el',num2str(1)] ;
    for j = 2:Kest
        str = ['el',num2str(j)];
        if length(str)>size(compNames,2)
            compNames = [compNames repmat(' ',j-1,1)];
        end
        compNames(j,:) = str;
    end
    misc(1).compNames=cellstr(compNames);
    
    %      misc.baseline = [];
    %      misc.residual =  filrecon-filoriginal;
    %      misc.filter_information.order = Kest;
    %      misc.individual_components = individual_components;

end
misc(1).begin = begin;
misc(1).ndp = ndp;
misc(1).frequency = frequency;
misc(1).step = step;

procres.misc = misc;
procres.scores = scores;
savef = get(handles.editsaveto,'String');
save(savef,'procres','step','begin','frequency');
% save(savef,'filreconAll','baselineAll','filoriginalAll','filter_information',...
%             'amplAll','dampAll','phasAll','freqAll',...
%             'individual_componentsAll','residualAll','compNames','frequency','step','begin');



% --- Executes on selection change in method.
function method_Callback(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method


% --- Executes during object creation, after setting all properties.
function method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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


function editorder_Callback(hObject, eventdata, handles)
% hObject    handle to editorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editorder as text
%        str2double(get(hObject,'String')) returns contents of editorder as a double


% --- Executes during object creation, after setting all properties.
function editorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editorder (see GCBO)
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

