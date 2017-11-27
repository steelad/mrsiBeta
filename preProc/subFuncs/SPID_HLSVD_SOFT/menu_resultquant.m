function varargout = menu_resultquant(varargin)
% RESULTQUANT M-file for resultquant.fig
%      RESULTQUANT, by itself, creates a new RESULTQUANT or raises the existing
%      singleton*.
%
%      H = menu_resultquant returns the handle to a new menu_resultquant or the handle to
%      the existing singleton*.
%
%      menu_resultquant('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in menu_resultquant.M with the given input arguments.
%
%      menu_resultquant('Property','Value',...) creates a new menu_resultquant or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before menu_resultquant_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to menu_resultquant_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help menu_resultquant

% Last Modified by GUIDE v2.5 28-Nov-2007 09:02:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @menu_resultquant_OpeningFcn, ...
    'gui_OutputFcn',  @menu_resultquant_OutputFcn, ...
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


% --- Executes just before menu_resultquant is made visible.
function menu_resultquant_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to menu_resultquant (see VARARGIN)

% Choose default command line output for menu_resultquant
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes menu_resultquant wait for user response (see UIRESUME)
% uiwait(handles.menuplot);
if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = menu_resultquant_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function initialize_gui(fig_handle, handles)
loadresults(fig_handle, handles);

function loadresults(fig_handle, handles)
cdir = getappdata(0,'currentdir');
[fn p]=uigetfile([cdir,'*.mat'],'Load result file from');
if isequal(fn,0) | isequal(p,0)
    return
end
cd(pwd);
loadf = fullfile(p,fn);
setappdata(0, 'currentdir',p);
setappdata(0, 'quantresultfile',loadf)

update(fig_handle, [], handles,loadf);
% if ~isfield(astr,'amplAll')|~isfield(astr,'compNames')
%     errordlg('The variables ''amplAll'' and ''compNames'' are not present in the mat file.');
% end





function update(hObject, eventdata, handles,loadf)
astr = load(loadf);

if isfield(astr,'procres')
    if ~isfield(astr.procres,'misc')
        if ~isfield(astr.procres,'filreconAll')|~isfield(astr.procres,'baselineAll')|...
                ~isfield(astr.procres,'step')|~isfield(astr.procres,'begin')|...
                ~isfield(astr.procres,'amplAll')|~isfield(astr.procres,'freqAll')|~isfield(astr.procres,'dampAll')|...
                ~isfield(astr.procres,'phasAll')|~isfield(astr.procres,'compNames')
            errordlg('Not enough variables or variables misspelled. Please check.');
            return
        end
    else
        miscOK = 1;
    end
else miscOK = 0;
end


if ~isfield(astr,'DBnorm')
    DBnorm = 0;
end


if miscOK
    for i = 1 :length(astr.procres.misc)
        amplAll(i,:) = astr.procres.misc(i).ampl;
        freqAll(i,:) = astr.procres.misc(i).freq;
        phasAll(i,:) = astr.procres.misc(i).phas;
        dampAll(i,:)= astr.procres.misc(i).damp;
    end
    i = 1;
    begin = astr.procres.misc(i).begin;
    step = astr.procres.misc(i).step;
    frequency = astr.procres.misc(i).frequency;
    ndp = astr.procres.misc(i).ndp;
    compNames = astr.procres.misc(i).compNames;
    names = fieldnames(astr.procres.misc);
    m=0;
    for ii=1:length(names)
        v = genvarname(names{ii});
        if ~strcmp(v,'compNames')
            for j = 1:length(astr.procres.misc)
                mt = size(astr.procres.misc(j).(names{ii}),2);
                if mt>m
                    m = mt;
                end
            end
            for j = 1:length(astr.procres.misc)
                b = astr.procres.misc(j).(names{ii});
                if ~isa(b,'struct')
                    if j==1
                        a = size(b,2);
                        if a <m
                            s =[b,zeros(size(b,1),m-a)];
                            eval([v ' = s;']);
                        else
                            eval([v ' = b(:,1:m);']);
                        end
                    else
                        a = size(b,2);
                        if a <m
                            s =[b,zeros(size(b,1),m-a)];
                            eval([v ' = [' v '; s];']);
                        else
                            eval([v ' = [' v '; b(:,1:m)];']);
                        end
                    end
                end
            end
        end
        m = 0;
    end
    ss=size(filoriginal,1);
    ndp = size(filoriginal,2);
    if ss==0
        %    set(handles.plotbutton,'Enable','off');
    elseif ss<=1
        set(handles.slider1,'Enable','off');
    else
        set(handles.slider1,'Enable','on');
        %    set(handles.plotbutton,'Enable','on');
        set(handles.slider1,'Min',1,'Max',ss,'SliderStep',[1/(ss-1) 2/(ss-1)]);
        nbr = str2double(get(handles.editspectrum,'string'));
        set(handles.slider1, 'Value',nbr);
        if get(handles.slider1,'Value')>ss
            set(handles.slider1, 'Value', 1);
        end
    end
else
    names = fieldnames(astr);
    begin = astr.begin;
    step = astr.step;

    for i=1:length(names)
        v = genvarname(names{i});
        eval([v ' = astr.(names{i});']);
    end
    ss=size(filoriginalAll,1);
    ndp = size(filoriginalAll,2);
    if ss==0
        %    set(handles.plotbutton,'Enable','off');
    elseif ss<=1
        set(handles.slider1,'Enable','off');
    else
        set(handles.slider1,'Enable','on');
        %    set(handles.plotbutton,'Enable','on');
        set(handles.slider1,'Min',1,'Max',ss,'SliderStep',[1/(ss-1) 2/(ss-1)]);
        nbr = str2double(get(handles.editspectrum,'string'));
        set(handles.slider1, 'Value',nbr);
        if get(handles.slider1,'Value')>ss
            set(handles.slider1, 'Value', 1);
        end
    end

    if ischar(compNames)
        for i = 1: size(compNames,1)% we suppose compNames a char array s.t. each row corresponds to one metabolite name
            t{i} = compNames(i,:);
        end
        compNames = t; % compNames is converted from char array to cell
    end

    if ~DBnorm

        if length(compNames) ~= size(amplAll,2)
            amplAll = amplAll';
            if length(compNames) ~= size(amplAll,2)
                errordlg('The number of metabolites of ''compNames'' and ''amplAll'' does not match.')
            end
        end
        if length(compNames) ~= size(dampAll,2)
            dampAll = dampAll';
            if length(compNames) ~= size(dampAll,2)
                errordlg('The number of metabolites of ''compNames'' and ''dampAll'' does not match.')
            end
        end
        if length(compNames) ~= size(freqAll,2)
            freqAll = freqAll';
            if length(compNames) ~= size(freqAll,2)
                errordlg('The number of metabolites of ''compNames'' and ''freqAll'' does not match.')
            end
        end
        if length(compNames) ~= size(phasAll,2)
            phasAll = phasAll';
            if length(compNames) ~= size(phasAll,2)
                errordlg('The number of metabolites of ''compNames'' and ''phasAll'' does not match.')
            end
        end
    end
end
nbr = str2double(get(handles.editspectrum,'string'));
% new text label parameters
xoffset = 0.02; yoffset = 0.9; vspedit = 0.014; c=0; textwidth = 0.09; textheight = 0.014;
hoffset = 0.11;
bg = [1 1 1];
nt  =size(amplAll,2);
i=1;
n = ['text', num2str(i)];
% eval(['o =' handles.(n)])
% while exist('[''handles.text'', num2str(i)]')
handlesplot = guihandles(gcf);
while isfield(handlesplot,n)
    %     cla(handlesplot.(n));
    set(handlesplot.(n),'string','');
    i = i+1;
    n = ['text', num2str(i)];
    %     eval(['o =' handles.n;])
end

handles.text101 = uicontrol('style','text','Units','normalized',...
    'Position',[xoffset+hoffset yoffset-c*vspedit textwidth textheight],...
    'HorizontalAlignment','center','String','Amplitude (a.u.)','BackgroundColor',bg);
handles.text102 = uicontrol('style','text','Units','normalized',...
    'Position',[xoffset+2*hoffset yoffset-c*vspedit textwidth textheight],...
    'HorizontalAlignment','center','String','Damping (kHz.rad)','BackgroundColor',bg);
c=c+1;

if DBnorm
    str1 = ['text' num2str(((i-1)*3)+1)];
    str2 = ['text' num2str(((i-1)*3)+2)];
    str3 = ['text' num2str(((i-1)*3)+3)];
    handles.(str1) = uicontrol('style','text','Units','normalized',...
        'Position',[xoffset yoffset-c*vspedit textwidth textheight],...
        'HorizontalAlignment','center','String',char(compNames{nbr}),'Tag',str1,'BackgroundColor',bg);
    handles.(str2) = uicontrol('style','text','Units','normalized',...
        'Position',[xoffset+hoffset yoffset-c*vspedit textwidth textheight],...
        'HorizontalAlignment','center','String',amplAll(nbr,1),'Tag',str2);
    %     set(handles.(str),'String',num2str(amplAll(nbr,i)));

    handles.(str3) = uicontrol('style','text','Units','normalized',...
        'Position',[xoffset+2*hoffset yoffset-c*vspedit textwidth textheight],...
        'HorizontalAlignment','center','String',dampAll(nbr,1),'Tag',str3);
    %     set(handles.(str),'String',num2str(dampAll(nbr,i)));
    c=c+1;
else
    for i = 1 : nt
        str1 = ['text' num2str(((i-1)*3)+1)];
        str2 = ['text' num2str(((i-1)*3)+2)];
        str3 = ['text' num2str(((i-1)*3)+3)];
        handles.(str1) = uicontrol('style','text','Units','normalized',...
            'Position',[xoffset yoffset-c*vspedit textwidth textheight],...
            'HorizontalAlignment','center','String',char(compNames{i}),'Tag',str1,'BackgroundColor',bg);
        handles.(str2) = uicontrol('style','text','Units','normalized',...
            'Position',[xoffset+hoffset yoffset-c*vspedit textwidth textheight],...
            'HorizontalAlignment','center','String',amplAll(nbr,i),'Tag',str2);
        %     set(handles.(str),'String',num2str(amplAll(nbr,i)));

        handles.(str3) = uicontrol('style','text','Units','normalized',...
            'Position',[xoffset+2*hoffset yoffset-c*vspedit textwidth textheight],...
            'HorizontalAlignment','center','String',dampAll(nbr,i),'Tag',str3);
        %     set(handles.(str),'String',num2str(dampAll(nbr,i)));
        c=c+1;
    end
end

c=c+1;
handles.text103 = uicontrol('style','text','Units','normalized',...
    'Position',[xoffset+hoffset yoffset-c*vspedit textwidth textheight],...
    'HorizontalAlignment','center','String','Frequency (kHz)','BackgroundColor',bg);
handles.text104 = uicontrol('style','text','Units','normalized',...
    'Position',[xoffset+2*hoffset yoffset-c*vspedit textwidth textheight],...
    'HorizontalAlignment','center','String','Phase (deg)','BackgroundColor',bg);
c=c+1;
if DBnorm
    str4 = ['text' num2str(((nt+i-1)*3)+1)];
    str5 = ['text' num2str(((nt+i-1)*3)+2)];
    str6 = ['text' num2str(((nt+i-1)*3)+3)];
    handles.(str4) = uicontrol('style','text','Units','normalized',...
        'Position',[xoffset yoffset-c*vspedit textwidth textheight],...
        'HorizontalAlignment','center','String',char(compNames{nbr}),'Tag',str4,'BackgroundColor',bg);

    handles.(str5) = uicontrol('style','text','Units','normalized',...
        'Position',[xoffset+hoffset yoffset-c*vspedit textwidth textheight],...
        'HorizontalAlignment','center','String',freqAll(nbr,1),'Tag',str5);

    handles.(str6) = uicontrol('style','text','Units','normalized',...
        'Position',[xoffset+2*hoffset yoffset-c*vspedit textwidth textheight],...
        'HorizontalAlignment','center','String',phasAll(nbr,1),'Tag',str6);
    %'HorizontalAlignment','center','String',360/2/pi*phasAll(nbr,1),'Tag',str6);
    c=c+1;
else
    for i = 1:nt
        str4 = ['text' num2str(((nt+i-1)*3)+1)];
        str5 = ['text' num2str(((nt+i-1)*3)+2)];
        str6 = ['text' num2str(((nt+i-1)*3)+3)];
        handles.(str4) = uicontrol('style','text','Units','normalized',...
            'Position',[xoffset yoffset-c*vspedit textwidth textheight],...
            'HorizontalAlignment','center','String',char(compNames{i}),'Tag',str4,'BackgroundColor',bg);

        handles.(str5) = uicontrol('style','text','Units','normalized',...
            'Position',[xoffset+hoffset yoffset-c*vspedit textwidth textheight],...
            'HorizontalAlignment','center','String',freqAll(nbr,i),'Tag',str5);

        handles.(str6) = uicontrol('style','text','Units','normalized',...
            'Position',[xoffset+2*hoffset yoffset-c*vspedit textwidth textheight],...
            'HorizontalAlignment','center','String',phasAll(nbr,i),'Tag',str6);
                  %  'HorizontalAlignment','center','String',360/2/pi*phasAll(nbr,i),'Tag',str6);
        c=c+1;
    end
end
plotinaxis(hObject, eventdata, handles,loadf)




function plotinaxis(hObject, eventdata, handles,loadf)
astr = load(loadf);
if isfield(astr,'procres')
    if isfield(astr.procres,'misc')
        miscOK = 1;
        begin = astr.procres.misc(1).begin;
        step = astr.procres.misc(1).step;
        frequency = astr.procres.misc(1).frequency;
        ndp = astr.procres.misc(1).ndp;
        compNames = astr.procres.misc(1).compNames;
        names = fieldnames(astr.procres.misc);
        m = 0;
        for ii=1:length(names)
            v = genvarname(names{ii});
            if ~strcmp(v,'compNames')&~strcmp(v,'step')&~strcmp(v,'begin')&~strcmp(v,'frequency')&~strcmp(v,'ndp')%&~strcmp(v,'individual_components')
                for j = 1:length(astr.procres.misc)
                    mt = size(astr.procres.misc(j).(names{ii}),2);
                    if mt>m
                        m = mt;
                    end
                end
                for j = 1:length(astr.procres.misc)
                    b = astr.procres.misc(j).(names{ii});
                    if ~isa(b,'struct')
                        if j==1
                            a = size(b,2);
                            if a <m
                                s =[b,zeros(size(b,1),m-a)];
                                eval([v ' = s;']);
                            else
                                eval([v ' = b(:,1:m);']);
                            end
                        else
                            a = size(b,2);
                            if a <m
                                s =[b,zeros(size(b,1),m-a)];
                                eval([v ' = [' v '; s];']);
                            else
                                eval([v ' = [' v '; b(:,1:m)];']);
                            end
                        end
                    end
                end
            end
            m = 0;
        end
    else
        begin = astr.begin;
        step = astr.step;
        names = fieldnames(astr);
        for i=1:length(names)
            v = genvarname(names{i});
            eval([v ' = astr.(names{i});'])
        end
        miscOK = 0;
    end
    % reshaping ----------------------------------
    if exist('filoriginalAll')
        if size(filoriginalAll,1)>size(filoriginalAll,2)
            filoriginalAll = filoriginalAll';
        end
    end
    if exist('filreconAll')
        if size(filreconAll,1)>size(filreconAll,2)
            filreconAll = filreconAll';
        end
    end
    if exist('baselineAll')
        if size(baselineAll,1)>size(baselineAll,2)
            baselineAll = baselineAll';
        end
    end
end
% ---------------------------------------------
if miscOK
    ndp = size(filoriginal,2);
    nbr = str2double(get(handles.editspectrum,'string'));
    fig = openfig('menu_resultquant', 'reuse');
    param.ndp = ndp;
    if exist('frequency')
        param.frequency = frequency;
    else param.frequency = 63856; % in kHz
        warning('The spectrometer frequency has been set at 63856 kHz.')
    end
    param.step = astr.procres.misc(1).step;
    param.begin = astr.procres.misc(1).begin;

    %plotting in top plot (residuals)
    axes(handles.axes1);
    if exist('filoriginal')&exist('filrecon')%&exist('baselineAll')
        filerr = filrecon(:,1:ndp)-filoriginal(:,1:ndp);
        displayplot(fig,[],guidata(fig),filerr,'nokeep',0,1,0,param);
    end
    title('residuals (filtered estimated signal - filtered original signal)')

    %plotting in second plot (baseline)
    axes(handles.axes2);
    if exist('baseline')
        if size(baseline,2)~=0
            displayplot(fig,[],guidata(fig),baseline(:,1:ndp),'nokeep',0,1,0,param);
        end
    end
    title('baseline')

    % %plotting in third plot (components)
    axes(handles.axes3);
    dimension = 3; %plot in 3D
    if exist('individual_components')
        param.ndp = size(individual_components,2);
        displayplot(fig,[],guidata(fig),individual_components((nbr-1)*size(ampl,2)+1:nbr*size(ampl,2),:),'nokeep',0,1,0,param,dimension);
    end
    title('individual corrected components')

    %plotting in bottom plot (filtered original vs filtered reconstructed estimated signal)
    axes(handles.axes4);
    param.ndp = ndp;
    if exist('filoriginal')&exist('filrecon')%&exist('baselineAll')
        displayplot(fig,[],guidata(fig),filoriginal(:,1:ndp),'nokeep',0,1,0,param);
        displayplot(fig,[],guidata(fig),filrecon(:,1:ndp),'keep',0,1,0,param);
    end
    title('filtered estimated signal  - filtered original signal (blue)')
end


% ---------------------------------------------
if ~miscOK
    ndp = size(filoriginalAll,2);
    nbr = str2double(get(handles.editspectrum,'string'));
    fig = openfig('menu_resultquant', 'reuse');
    param.ndp = ndp;
    if exist('frequency')
        param.frequency = frequency;
    else param.frequency = 63856; % in kHz
        warning('The spectrometer frequency has been set at 63856 kHz.')
    end
    param.step = step;
    param.begin = begin;

    %plotting in top plot (residuals)
    axes(handles.axes1);
    if exist('filoriginalAll')&exist('filreconAll')%&exist('baselineAll')
        % filerr = filreconAll(nbr,1:ndp)-filoriginalAll(nbr,1:ndp);
        %filerr = filreconAll(:,1:ndp)+baselineAll(:,1:ndp)-filoriginalAll(:,1:ndp);
        filerr = filreconAll(:,1:ndp)-filoriginalAll(:,1:ndp);
        displayplot(fig,[],guidata(fig),filerr,'nokeep',0,1,0,param);
    end
    title('residuals (filtered estimated signal - filtered original signal)')

    %plotting in second plot (baseline)
    axes(handles.axes2);
    if exist('baselineAll')
        if size(baselineAll,2)~=0
            % displayplot(fig,[],guidata(fig),baselineAll(nbr,1:ndp),'nokeep',0,1,0,param);
            displayplot(fig,[],guidata(fig),baselineAll(:,1:ndp),'nokeep',0,1,0,param);
        end
    end
    title('baseline')

    % %plotting in third plot (components)
    axes(handles.axes3);
    dimension = 3; %plot in 3D
    if exist('individual_componentsAll')
        param.ndp = size(individual_componentsAll(nbr).indcomp,2);
        % displayplot(fig,[],guidata(fig),individual_componentsAll(nbr).indcomp(:,1:ndp),'nokeep',0,1,0,param,dimension);
        displayplot(fig,[],guidata(fig),individual_componentsAll(nbr).indcomp,'nokeep',0,1,0,param,dimension);
    end
    title('individual corrected components')

    %plotting in bottom plot (filtered original vs filtered reconstructed estimated signal)
    axes(handles.axes4);
    param.ndp = ndp;
    if exist('filoriginalAll')&exist('filreconAll')%&exist('baselineAll')
        displayplot(fig,[],guidata(fig),filoriginalAll(:,1:ndp),'nokeep',0,1,0,param);
        displayplot(fig,[],guidata(fig),filreconAll(:,1:ndp),'keep',0,1,0,param);
        %displayplot(fig,[],guidata(fig),baselineAll(:,1:ndp),'keep',0,1,0,param);
        % displayplot(fig,[],guidata(fig),filoriginalAll(nbr,1:ndp),'nokeep',0,1,0,param);
        % displayplot(fig,[],guidata(fig),filreconAll(nbr,1:ndp),'keep',0,1,0,param);
    end
    title('filtered estimated signal  - filtered original signal (blue)')
end


function plotting2D(fig_handle, handles)

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.editspectrum,'String',num2str(round(get(handles.slider1,'value'))));
%plotinaxis(hObject, eventdata, handles,getappdata(0, 'quantresultfile'));
update(hObject, eventdata, handles,getappdata(0, 'quantresultfile'))

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in plotfreq.
function plotfreq_Callback(hObject, eventdata, handles)
% hObject    handle to plotfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plotfreq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotfreq
update(hObject, eventdata, handles,getappdata(0, 'quantresultfile'));

% --- Executes during object creation, after setting all properties.
function plotfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in plotreal.
function plotreal_Callback(hObject, eventdata, handles)
% hObject    handle to plotreal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plotreal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotreal
update(hObject, eventdata, handles,getappdata(0, 'quantresultfile'));

% --- Executes during object creation, after setting all properties.
function plotreal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotreal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editspectrum_Callback(hObject, eventdata, handles)
% hObject    handle to editspectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editspectrum as text
%        str2double(get(hObject,'String')) returns contents of editspectrum as a double
set(handles.slider1,'value',round(str2double(get(handles.editspectrum,'string'))));
update(hObject, eventdata, handles,getappdata(0, 'quantresultfile'));

% --- Executes during object creation, after setting all properties.
function editspectrum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editspectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --------------------------------------------------------------------
function loadresults_Callback(hObject, eventdata, handles)
% hObject    handle to loadresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadresults(hObject, handles);

% --- Executes during object creation, after setting all properties.
function spectext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --------------------------------------------------------------------
function savefig_Callback(hObject, eventdata, handles)
% hObject    handle to savefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savefigure(hObject, eventdata, handles);



