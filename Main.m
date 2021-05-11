%%% -----------------------------------------------------------------------
clc, clearvars, close all;
%%% -----------------------------------------------------------------------
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = ((0:L-1)*T)';     % Time vector
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
X(:,1) = S + .2*randn(size(t));
X(:,2) = S + .3*randn(size(t));
%%% -----------------------------------------------------------------------
SigPro_Input.ti    = [];
SigPro_Input.yi    = X;
SigPro_Input.fs    = Fs;
SigPro_Input.Wndw  = [];
SigPro_Input.Trend = 0;
SigPro_Input.ts    = [];
SigPro_Input.ffi   = 0.1;
SigPro_Input.fff   = 20;
SigPro_Input.plotts  = 'Acel. (m/s^2)';
SigPro_Input.plotfft = 0;
SigPro_Input.plotpsd = 0;
SigPro_Input.plotspt = 0;
%%% -----------------------------------------------------------------------
SignOut = SigPro(SigPro_Input.ti,...
    SigPro_Input.yi,...
    SigPro_Input.fs,...
    SigPro_Input.Wndw,...
    SigPro_Input.Trend,...
    SigPro_Input.ts,...
    SigPro_Input.ffi,...
    SigPro_Input.fff,...
    SigPro_Input.plotts,...
    SigPro_Input.plotfft,...
    SigPro_Input.plotpsd,...
    SigPro_Input.plotspt);
%%% -----------------------------------------------------------------------