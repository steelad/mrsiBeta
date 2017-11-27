function varargout = menu_DBnorm(varargin)
% MENU_DBNORM M-file for menu_DBnorm.fig
%      MENU_DBNORM, by itself, creates a new MENU_DBNORM or raises the existing
%      singleton*.
%
%      H = MENU_DBNORM returns the handle to a new MENU_DBNORM or the handle to
%      the existing singleton*.
%
%      MENU_DBNORM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MENU_DBNORM.M with the given input arguments.
%
%      MENU_DBNORM('Property','Value',...) creates a new MENU_DBNORM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before menu_DBnorm_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to menu_DBnorm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help menu_DBnorm

% Last Modified by GUIDE v2.5 28-Feb-2007 09:35:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @menu_DBnorm_OpeningFcn, ...
    'gui_OutputFcn',  @menu_DBnorm_OutputFcn, ...
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
set(handles.browserf2, 'UserData', handles.editloadfrom2);
set(handles.browserp, 'UserData', handles.editsaveto);
set(handles.editfilterlength,'String',50);


% --- Executes just before menu_DBnorm is made visible.
function menu_DBnorm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to menu_DBnorm (see VARARGIN)

% Choose default command line output for menu_DBnorm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes menu_DBnorm wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = menu_DBnorm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in normalize.
function normalize_Callback(hObject, eventdata, handles)
% hObject    handle to normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

typemethode = get(handles.quantmethod,'Value'); %typemethode = 1 if AQSES
%         = 2 if HLSVD-PRO

filtertype = get(handles.editfiltertype,'Value'); %  filtertype = 1 if MP-FIR
%        = 2 if HLSVD-PRO

ff  = get(handles.editloadfrom,'String');
if strcmp(ff,'')
    errordlg('You must enter a database file');
    return
else
    db = load(ff);
    if ~isfield(db,'begin')|~isfield(db,'step')|~isfield(db,'ndp')
        errordlg('The parameters begin, step and ndp do not exist. Please add them in the mat file of your database.')
        return
    end

    name_metabolite = get(handles.editmetabolites,'String');%{'Myo','Pch','Cr','Glu','NAA','Lac','Lip1','Lip2'};
    nbrm = str2num(get(handles.editnbrmet,'String'));
    ffread(name_metabolite);
    w = ffread(nbrm);
    %     figure
    %     plot(real(eval(['db.',char(w(3))])))
    if nbrm == 1
        dbase_landmark1(1,:)= eval(['db.',char(w)]);
    else
        for i=1:length(w)
            dbase_landmark1(i,:) = eval(['db.',char(w(i))]);
        end
    end
    dbase_landmark1 = dbase_landmark1'; % each metabolite profile is a column
end

ff  = get(handles.editloadfrom2,'String');
if strcmp(ff,'')
    errordlg('You must enter a database file');
    return
else
    db2 = load(ff);
    if ~isfield(db2,'begin')|~isfield(db2,'step')|~isfield(db2,'ndp')
        errordlg('The parameters begin, step and ndp do not exist. Please add them in the mat file of your database.')
        return
    end
    begin = db2.begin;
    step = db2.step;
    ndp = db2.ndp;
    frequency = db2.frequency;
    lambda = 0.01/step;
    set(handles.editallowdamp,'String',lambda*10);
    set(handles.editallowfreq,'String',lambda*10);
    t   = [0:step:(ndp-1)*step] + begin;%in ms
    if db.ndp~=ndp
        diff = ndp-db.ndp;
        if diff<0
            warning('The database has been truncated to fit the number of points of the signal')
            dbase_landmark1 = dbase_landmark1(:,1:ndp);
        else
            warning('The database has been extended by', diff,' to fit the number of points of the signal')
            dbase_landmark1 = [dbase_landmark1 zeros(size(dbase_landmark1,1),diff)];
        end
    end
end


% Filtering
if  filtertype==1
    filtertype = '';
elseif filtertype==2
    filtertype = 'fir';
elseif filtertype==3
    filtertype = 'hlsvd';
end

allow_damp = str2num(get(handles.editallowdamp,'String'));
allow_freq = str2num(get(handles.editallowfreq,'String'))*2*pi;
lambda = str2num(get(handles.editlambda,'String'));
maxiter = str2num(get(handles.editmaxiter,'String'));
phasedistort = str2num(get(handles.editECC,'String'));


ppmrange = str2num(get(handles.editrange,'String'));
fsampl = 1/step;
xmin_mouse = ((ppmrange(1)-4.7)/(fsampl/(frequency/1000)))/1000;%in kHz
xmax_mouse = ((ppmrange(2)-4.7)/(fsampl/(frequency/1000)))/1000;%in kHz
ripple = str2num(get(handles.editripple,'String'));
filterlength = str2num(get(handles.editfilterlength,'String'));
equalph = get(handles.checkequalph,'Value'); 
baseline_boolean = get(handles.checkbaseline,'Value');

if get(handles.checkplotiter,'Value')
    plotting = 1;
else
    plotting = 0;
end


dbase_landmark1temp = dbase_landmark1;% each metabolite profile is a column
for i = 1:length(w)
    dbase_landmark1 = dbase_landmark1temp(:,i);
    if nbrm == 1
        signal = eval(['db2.',char(w)]);
    else
        signal = eval(['db2.',char(w(i))]);
    end
    abs_scl_factor = ones(size(dbase_landmark1,2),1);
    index1 = [1:size(dbase_landmark1,2)];
    nrstofs = size(dbase_landmark1,2);
    squares = zeros(nrstofs',1);

    filtertmp = filtertype;
    name_metplot = get(handles.editmetplot,'String');
    [A, count, errmsg, nextindex] = sscanf(get(handles.editmetplot,'String'),'%s');
    nbrm2 = count;
    if count
        ffread(name_metplot);
        w2 = ffread(nbrm2);
        if count==1 & nbrm==1
            if strcmp(w2,char(w))
                filtertype = '';
            end
        elseif count==1 & nbrm~=1
            if strcmp(w2,char(w(i)))
                filtertype = '';
            end
        elseif count~=1 & nbrm==1
            for ii = 1:length(w2)
                if strcmp(char(w),char(w2(ii)))
                    filtertype = '';
                    break
                end
            end
        elseif count~=1 & nbrm~=1

            for ii = 1:length(w2)
                if strcmp(char(w(i)),char(w2(ii)))
                    filtertmp = filtertype;
                    filtertype = '';
                    break
                end
            end
        end
    end

    %prior knowledge
    prior_equal = '';
    equal_to = '';


    % initial values for nonlinear parameters
    zer    = zeros(size(index1));
    demp0  = zer;
    gauss0 = zer;
    freq0  = zer;
    e0     = zer;

    j=sqrt(-1);
    dbase_landmark1 = real(dbase_landmark1)-j*(imag(dbase_landmark1));
    % figure
    % plot(real(fftshift(fft(dbase_landmark1(:,3)))))
    % figure
    % plot(real(dbase_landmark1(:,3)))

    if (lambda<0)

        [recon,baseline,filrecon,filbaseline,filoriginal,filter_information,...
            ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,CR_error,...
            individual_components,residual] = ...
            akses_gcv(signal,dbase_landmark1,t,step,abs_scl_factor,...
            demp0,freq0,gauss0,e0,1,0,0,phasedistort,...
            filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
            prior_equal,equal_to,allow_damp,allow_freq,baseline_boolean,lambda,...
            maxiter,plotting);
    else

        [recon,baseline,filrecon,filbaseline,filoriginal,filter_information,...
            ampl,demp,phas,freq,gauss,e,ampl_ratios,ampl_absolute,...
            CR_error,individual_components,residual] = ...
            akses_fit(signal,dbase_landmark1,t,step,abs_scl_factor,...
            demp0,freq0,gauss0,e0,1,0,0,phasedistort,...
            filtertype,xmin_mouse,xmax_mouse,ripple,filterlength,...
            prior_equal,equal_to,allow_damp,allow_freq,equalph,baseline_boolean,lambda,...
            maxiter,plotting);

    end
    

        reconAll(i,:) = (real(recon)-sqrt(-1)*imag(recon))';
        baselineAll(i,:) = (real(baseline)-sqrt(-1)*imag(baseline))';
        if i>1
        if length(filrecon)>size(filreconAll,2)
            filreconAll = [filreconAll,zeros(size(filreconAll,1),length(filrecon)-size(filreconAll,2))];
            filbaselineAll = [filbaselineAll,zeros(size(filbaselineAll,1),length(filrecon)-size(filreconAll,2))];
            filoriginalAll = [filoriginalAll,zeros(size(filoriginalAll,1),length(filrecon)-size(filreconAll,2))];
%             individual_components(i).individual_components = [individual_components, zeros(size(signal,1),size(dbase_landmark1,2),length(filrecon)-size(filreconAll,2))];
        end
        if length(residual)>size(residualAll,2)
            residualAll =  [residualAll,zeros(size(residualAll,1),length(residual)-size(residualAll,2))];
        end
        end
        filreconAll(i,1:length(filrecon)) = (real(filrecon)-sqrt(-1)*imag(filrecon))';
        filbaselineAll(i,1:length(filrecon)) = filbaseline';
        filoriginalAll(i,1:length(filrecon)) = (real(filoriginal)-sqrt(-1)*imag(filoriginal))';
        residualAll(i,1:length(residual)) = (real(residual)-sqrt(-1)*imag(residual))';
%     end
    
    amplAll(i,:) = ampl';
    dampAll(i,:) = demp';
    phasAll(i,:) = phas';
    freqAll(i,:) = freq/(2*pi);
    gaussAll(i,:) = gauss;
    eAll(i,:) = e';
    ampl_ratiosAll(i,:)=ampl_ratios';
    ampl_absoluteAll(i,:) = ampl_absolute';
    CRname = ['CR_error',num2str(i)];
    CRAll.(CRname)  = CR_error;
    
    %         ic = individual_components;
    %         clear individual_components
    individual_componentsAll(i).indcomp = individual_components;
    filtertype = filtertmp;
end %for i = 1:length(w)

amplsimul = getappdata(0,'amplsimul');
freqsimul = getappdata(0,'freqsimul');
phassimul = getappdata(0,'phassimul');
dampsimul = getappdata(0,'dampsimul');

if strcmp(amplsimul,'')|strcmp(dampsimul,'')|strcmp(freqsimul,'')|strcmp(phassimul,'')
    disp('WARNING: The data do not contain information about the simulated parameters. You will not be able to check the exact accuracy of your estimates.')
end

compNames = w;
DBnorm = 1;
savef = get(handles.editsaveto,'String');
if strcmp(savef,'')
    errordlg('You must first specify a file in which you want to store the results.')
else
    savef2 = [savef(1:end-4),'_res']; 
    save(savef2,'DBnorm','reconAll','baselineAll','filreconAll','filbaselineAll','filoriginalAll','filter_information',...
        'amplAll','dampAll','phasAll','freqAll','gaussAll','eAll','ampl_ratiosAll','ampl_absoluteAll','CRAll',...
        'individual_componentsAll','residualAll','amplsimul','dampsimul','freqsimul','phassimul','compNames','frequency','step','begin');
end

%Database normalization
save(savef,'begin','frequency','step','ndp')
for i = 1:length(w)
    if nbrm == 1
        v = genvarname(w);
        signal = db.(w);
        signal = signal*amplAll(i);
        eval([v ' = signal;'])
        save(savef,v,'-append')
    else
        v = genvarname(w{i});
        signal = db.(w{i});
        signal = signal*amplAll(i);
        eval([v ' = signal;'])
        save(savef,v,'-append')
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
% --- Executes on button press in browserf.
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



function editloadfrom2_Callback(hObject, eventdata, handles)
% hObject    handle to editloadfrom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editloadfrom2 as text
%        str2double(get(hObject,'String')) returns contents of editloadfrom2 as a double


% --- Executes during object creation, after setting all properties.
function editloadfrom2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editloadfrom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in browserf2.
function browserf2_Callback(hObject, eventdata, handles)
% hObject    handle to browserf2 (see GCBO)
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
%loadmetabolites(hObject, eventdata, handles);



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




% --- Executes on button press in checkbaseline.
function checkbaseline_Callback(hObject, eventdata, handles)
% hObject    handle to checkbaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbaseline


% --- Executes on button press in checkplotiter.
function checkplotiter_Callback(hObject, eventdata, handles)
% hObject    handle to checkplotiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkplotiter





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



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to editlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlambda as text
%        str2double(get(hObject,'String')) returns contents of editlambda as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
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



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to editmaxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmaxiter as text
%        str2double(get(hObject,'String')) returns contents of editmaxiter as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
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



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to editECC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editECC as text
%        str2double(get(hObject,'String')) returns contents of editECC as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
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


