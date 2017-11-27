function varargout = disp_concentrations(varargin)
% DISP_CONCENTRATIONS Application M-file for disp_concentrations.fig
%    FIG = DISP_CONCENTRATIONS launch disp_concentrations GUI.
%    DISP_CONCENTRATIONS('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 27-Jul-2004 13:38:43

if nargin == 0  % LAUNCH GUI

  fig = openfig(mfilename,'new');

  % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig);
  guidata(fig, handles);

  if nargout > 0
    varargout{1} = fig;
  end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
try
  
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
catch
  disp(lasterr);
end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


% --------------------------------------------------------------------
function varargout = disp_results_CreateFcn(h, eventdata, handles, recon,original, phas, step, t0, fir_h, xmin_mouse, xmax_mouse, baseline)
% Stub for CreateFcn of the axes handles.disp_results.

axes(h);
if (nargin<12) | (isempty(baseline)==1)
  baseline = zeros(length(recon),1);
end

ndft = length(baseline);

offset1 = disp_sigreal(baseline,mean(phas),step,t0,ndft,' ',0,3,fir_h); 
          hold on
offset2 = disp_sigreal(recon+baseline,mean(phas),step,t0,ndft,' ',...
                       offset1,2,fir_h);
offset3 = disp_sigreal(original,mean(phas),step,t0,ndft,' ',...
                       offset1,1,fir_h);
offset4 = max(offset2, offset3);
offset5 = disp_sigreal(recon,mean(phas),step,t0,ndft,' ',...
                       offset4,2,fir_h); 
offset6 = disp_sigreal(original-baseline,mean(phas),step,t0,ndft,' ',...
                       offset4,1,fir_h); 
disp_sigreal(original-recon-baseline,mean(phas),step,t0,ndft,' ',...
             max(offset5,offset6),1,fir_h); 

v = axis;
xmin_mouse=4.65+1000*xmin_mouse/63.86;
xmax_mouse=4.65+1000*xmax_mouse/63.86;
axis([xmin_mouse xmax_mouse v(3) v(4)]);  
grid, hold off  

%%%% disp_concentrations('disp_results_CreateFcn',gcbo,[],guidata(gcbo),recon,original, phas, step, t0, fir_h, xmin_mouse, xmax_mouse, baseline)

% --------------------------------------------------------------------
function varargout = Pushbutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Pushbutton1.
  close;



% --------------------------------------------------------------------
function varargout = print_button_Callback(h, eventdata, handles, title_analysis)
% Stub for Callback of the uicontrol handles.print_button.

  [FileName,PathName] = uiputfile('*.ps','Save ps file', 'test.ps');
  if isequal(FileName,0) | isequal(PathName,0)
       return
end
  if (FileName ~= 0)
    print('-dpsc2',  '-r600', '-adobecset', [PathName FileName])
  end
  

% --------------------------------------------------------------------
function varargout = show_individual_comp_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.show_individual_comp.

  show_individual_components
  %uiwait;


% --- Executes during object deletion, before destroying properties.
function disp_results_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to disp_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


