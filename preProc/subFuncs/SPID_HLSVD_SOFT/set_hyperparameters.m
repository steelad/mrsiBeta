function [load_saved_signal,baseline_boolean,allow_damp,allow_freq,phasedistort,...
    plotting,signal_type,bsline_type,noise_level,lambda,maxiter] = ...
    set_hyperparameters(argin);

% function that sets the flags and hyperparameter values


% Default parameters are: (you may edit them)

load_saved_signal = 1; % 1 to load a test example from file testdat/testnl,
% 2 to simulate random signal,
% other, to load another signal you might have
baseline_boolean  = 1; % 1 to fit a baseline, 0 otherwise
allow_damp      = .1;% variation allowed in damping and frequency shifts
allow_freq      = .1;% variation allowed in damping and frequency shifts
lambda            =.04;% value of regularization parameter
maxiter           = 25;% maximum number of iterations for nonlinear solver
phasedistort      = 0; % 1 - do Eddy current correction
% 0 - don't correct for Eddy current
plotting          = 1; % 1 for plotting at each iteration, 0 otherwise
signal_type   = 'ones';% 'ones' - simulate a signal such that all amplitudes
%         equal 1;
% 'rand' - simulate a signal with random amplitudes;
bsline_type = 'spline';% 'spline' - simulate baseline using splines;
% 'gauss'  - simulate baseline using several gaussians
noise_level    = 1e-12;% standard deviation of noise added to simulated signal


if (nargin > 0), return, end % stop if a the function is called
% with a dummy argument

% else, ask for user interaction

text = ['\n\nDefault parameters are:                    meaning:\n\n',...
    '1) load_saved_signal = 1;           %%   1 to load a test example from file testdat/testnl,\n',...
    '                                    %%   2 to simulate random signal,\n',...
    '                                    %%   other, to load another signal you might have\n',...
    '2) baseline_boolean  = 1;           %%   1 to fit a baseline, 0 otherwise\n',...
    '3) allow_damp      = .1;          %%   variation allowed in damping \n',...
    '4) phasedistort      = 0;           %%   1 - do Eddy current correction\n',...
    '                                    %%   0 - don''t correct for Eddy current\n',...
    '5) plotting          = 1;           %%   1 for plotting at each iteration, 0 otherwise\n',...
    '6) signal_type       = ''ones'';      %%   ''ones'' - simulate a signal such that all amplitudes equal 1;\n',...
    '                                    %%   ''rand'' - simulate a signal with random amplitudes;\n',...
    '7) bsline_type       = ''spline'';    %%   ''spline'' - simulate baseline using splines;\n',...
    '                                    %%   ''gauss''  - simulate baseline using several gaussians\n',...
    '8) noise_level       = 1e-12;       %%   standard deviation of noise added to simulated signal\n',...
    '9) lambda            =.04;          %%   value of regularization parameter  \n',...
    '10) maxiter          = 25;          %%   maximum number of iterations for nonlinear solver\n\n',...
'11) allow_freq      = .1;          %%   variation allowed in frequency shifts\n'];
    fprintf(text)

yes = input('Modify parameters? (y/n)  ','s');

if (length(yes)>0) & (strcmp(yes(1),'y') == 1)

    which = input('Which parameters? (Write their numbers between brackets [1,2,...])  ');

    if find(which == 1)
        load_saved_signal = input('load_saved_signal = (0,1,2) ');
    end
    if find(which == 2)
        baseline_boolean = input('baseline_boolean = (0,1) ');
    end
    if find(which == 3)
        allow_damp = input('allow_damp = (>0) ');
    end
    if find(which == 4)
        phasedistort = input('phasedistort = (0,1) ');
    end
    if find(which == 5)
        plotting = input('plotting = (0,1) ');
    end
    if find(which == 6)
        signal_type = input('signal_type = (ones, rand) ','s');
    end
    if find(which == 7)
        bsline_type = input('bsline_type = (spline, gauss) ','s');
    end
    if find(which == 8)
        noise_level = input('noise_level = (>0)');
    end

    if nargout > 8
        if find(which == 9)
            lambda = input('lambda = (>0)');
        end
        if find(which == 10)
            maxiter = input('maxiter = (int>0) ');
        end

        if find(which == 11)
            allow_freq = input('allow_freq = (>0) ');
        end
    end

end

