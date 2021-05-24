function SignOut = SigPro(varargin)
%
% SignOut = SigPro(ti,yi,fs,Wndw,trend,ts,ffi,fff,TimSigYLabel,plotfft,plotpsd,plotspt)
%
% Función para realizar el procesamiento de señales.
%
% INPUTS:
%
%   ti: puede ser [], escalar indicando en que tiempo inicial o
%       directamente el vector de tiempo inicial (antes de ser
%       remuestreado).
%
%   yi: señales iniciales (antes de ser procesadas y remuestreadas).
%
%   fs: frecuencia de muestreo, en caso de ser diferente de 0 y menor a la
%       del vector ti, se hará remuestreo.
%
%   Wndw: puede ser [] o vector con dos valores en *número de puntos* que limitan la ventana de
%       datos seleccionada. Ejemplo: Wndw = [1, 578]; % (seleccionar los
%       datos entre el segundo 1/fs y el 578/fs).
%
%   trend: valor lógico (1:true, 0:false) para realizar detrend.
%
%   ts: duración del suavizado inicial y final de la señal (eliminar ruido
%       cuando la excitación es nula).
%
%   ffi: frecuencia de corte del filtro pasa alto.
%
%   fff: frecuencia de corte del filtro pasa bajo.
%
%   TimSigYLabel: Character indicating YLabel for time-domain plot.
%
%   plotfft: 1:plot fft, 0:do not plot.
%
%   plotpsd: 1:plot psd, 0:do not plot.
%
%   plotspt: 1:plot spectrogram, 0:do not plot.
%
%   plotpwl: 1:plot pwelch, 0:do not plot.
%
% OUTPUTS:
%
%   SignOut: estructura con las propiedades de la señal procesada.
%
% Created By:
%
%   Mateo G. H. (2020/06/01)
%
% Modified By:
%
%   Mateo G. H. (2021/06/05)

Fcn_tic = tic;
%%% -----------------------------------------------------------------------
ti = varargin{1};
yi = varargin{2};
fs = varargin{3};
Wndw = varargin{4};
trend = varargin{5};
ts = varargin{6};
ffi = varargin{7};
fff = varargin{8};
TimSigYLabel = varargin{9};
plotfft = varargin{10};
plotpsd = varargin{11};
plotspt = varargin{12};
if nargin > 12
    plotpwl = varargin{13};
else
    plotpwl = 0;
end
%%%%%%%%%%%% --------------------------------------------------------------
%%% Wndw %%%
%%%%%%%%%%%%
if ~isempty(Wndw)
    yi = yi(Wndw(1):Wndw(end),:);
end
%%%%%%%%%%%%%%%% ----------------------------------------------------------
%%% Resample %%%
%%%%%%%%%%%%%%%%
if isempty(ti)
    fsi  = fs;
    tini = 0;
elseif length(ti)==1
    fsi  = fs;
    tini = ti;
else
    tini = 0;
    fsi  = round(1/(ti(2)-ti(1)),10);
end
if fs>0 && fs<fsi; y = resample(yi,fs,fsi); else; y = yi; fs = fsi; end
dt = 1/fs;
t  = round((0:(length(y(:,1))-1))'/fs+tini(1),10);
%%%%%%%%%%%%%%% -----------------------------------------------------------
%%% Detrend %%%
%%%%%%%%%%%%%%%
if trend==1
    y = detrend(y);
end
%%%%%%%%%%%%%% ------------------------------------------------------------
%%% Smooth %%%
%%%%%%%%%%%%%%
if ~isempty(ts)
    if ts>0
        S = ones(length(t),1);
        S(1:fs*ts/2-1) = (1-cos(2*pi*(1/ts)*t(1:fs*ts/2-1)))/2;
        S(end-(fs*ts/2-1):end) = (1+cos(2*pi*(1/ts)*t(1:fs*ts/2)))/2;
        S([1,end]) = [0;0];
        S = S*ones(size(y(1,:)));
        y = S.*y;
    end
end
%%%%%%%%%%%%%% ------------------------------------------------------------
%%% Filter %%%
%%%%%%%%%%%%%%
MaxFreq = min([0.4*fs,fff]);
if ffi~=0     	%%% High
    fc      = min([ffi,0.99*fs/2]); % (Hz)
    [z,p,k] = butter(6,fc/(fs/2),'high');
    sos     = zp2sos(z,p,k);
    y       = filtfilt(sos,1,y); % fvtool(sos,'Analysis','freq')
end
if fff~=0	%%% Low
    fc    = min([fff,0.99*fs/2]); % (Hz)
    [B,A] = butter(6,fc/(fs/2),'low');
    y     = filtfilt(B,A,y); % fvtool(B,A)
end
%%%%%%%%%%%% --------------------------------------------------------------
%%% Time %%%
%%%%%%%%%%%%
if ~strcmp(TimSigYLabel,' ') && ~isempty(TimSigYLabel)
    PlotTimeSignal(t,y,TimSigYLabel)
end
%%%%%%%%%%% ---------------------------------------------------------------
%%% FFT %%%
%%%%%%%%%%%
L  = length(y);             	% Longitud de la señal
Y  = 2*(fft(y)./L);             % (Y(f)) transformada rápida de Fourier.
Yp = angle(Y);                  % fase (rad).
f  = ((0:length(Y)/2-1)*fs/L)';	% (Hz)
Y  = Y(1:length(f),:);          % (|Y(f)|)
if plotfft==1
    PlotY(f,abs(Y),MaxFreq,'Frecuencia (Hz)','|Y(f)|','FFT');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  -------------------------------------------
%%% Power spectral density %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotpsd==1
    nfft  = length(y);
    windw = rectwin(nfft);
    [Ypsd,fpsd] = periodogram(y,windw,nfft,fs);
    PlotY(fpsd,10*log10(Ypsd),MaxFreq,'Frecuencia (Hz)','PSD (dB/Hz)','Periodogram power spectral density');
end
%%%%%%%%%%%%%%%%%%% -------------------------------------------------------
%%% Spectrogram %%%
%%%%%%%%%%%%%%%%%%%
if plotspt==1
    PlotSpectrogram(y,fs);
end
%%%%%%%%%%%%%% ------------------------------------------------------------
%%% PWelch %%%
%%%%%%%%%%%%%%
if plotpwl==1
    L = length(y);
    Nfft  = pow2(nextpow2(L));
    WndwL = round(L/20);
    Wndw  = hamming(WndwL);
    Nvrlp = floor(0.5*WndwL);
    [PW_Y,PW_f] = pwelch(y,Wndw,Nvrlp,Nfft,fs);
    %%% Plot
    PlotY(PW_f,10*log10(PW_Y),MaxFreq,'Frecuencia (Hz)','PSD (dB/Hz)','Welch''s power spectral density');
    %%%
    SignOut.Welch.f = PW_f;
    SignOut.Welch.Y = PW_Y;
end
%%%%%%%%%%%%%%%%%%% -------------------------------------------------------
%%% Export data %%%
%%%%%%%%%%%%%%%%%%%
SignOut.t   = t;
SignOut.y   = y;
SignOut.Y   = Y;
SignOut.f   = f;
SignOut.Yp  = Yp;
SignOut.fs  = fs;
SignOut.dt  = dt;
SignOut.ffi = ffi;
SignOut.fff = fff;
%%% -----------------------------------------------------------------------
disp(['SigPro: ',num2str(toc(Fcn_tic),'%.2f')])
end
%% PlotTimeSignal
function PlotTimeSignal(ti,yi,YLabel)
%
% PlotSignalCateg(Data,Case)
%
% INPUTS:
%   Data: structure obtained from SignalCateg plot.
%   Case: analysis window length as number of points.
%
% OUTPUTS:
%
% %%%%%%%%%%%%%%%%%%
% %%% Mateo G.H. %%%
% %%% 2021/03/15 %%%
% %%%%%%%%%%%%%%%%%%
GraphProp = GraphicalProperties;
%%% -----------------------------------------------------------------------
t = seconds(ti); t.Format = 'mm:ss';
A = yi;
%%% ---
wid  = 15;
hei  = 4;
Fig  = figure('Units','centimeters','Position',[1 6 wid hei]);
N    = size(A,2);
Clrs = colormap(winter(N));
LnSt = {'-','--','-.','-','--','-.'}; if N/length(LnSt)>1; LnSt = repmat(LnSt,1,ceil(N/length(LnSt))); end
%%% -----------------------------------------------------------------------
Ax1 = axes(Fig);
Pl1 = plot(t,A,'Color','k');
for ii=1:N
    Pl1(ii).Color     = Clrs(ii,:);
    Pl1(ii).LineWidth = GraphProp.linewidth;
    Pl1(ii).LineStyle = LnSt{ii};
end
%%% -----------------------------------------------------------------------
Ax1.XLim = [t(1),t(end)];
Ax1.YLim = [-1.1,1.1]*max(abs(A(:)));
Ax1.XLabel.String = 'Tiempo (mm:ss)';
Ax1.YLabel.String = YLabel;
set(Ax1,GraphProp.Prop);set(Ax1.XAxis,GraphProp.PropXA);set(Ax1.YAxis,GraphProp.PropYA);set(Ax1.Title,GraphProp.PropT);
set(Ax1.XLabel,GraphProp.PropXL);set(Ax1.YLabel,GraphProp.PropYL);
%%% -----------------------------------------------------------------------
lgn1 = cellstr(num2str((1:N)'));
lgn1 = legend(lgn1,'FontSize',GraphProp.fontsize-2,'Location','best');
set(lgn1.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.9]));
end
%%% -----------------------------------------------------------------------
%% PlotY
function PlotY(varargin)
GraphProp = GraphicalProperties;
%%% -----------------------------------------------------------------------
f = varargin{1};
Y = varargin{2};
MaxFreq = varargin{3};
XLabel = varargin{4};
YLabel = varargin{5};
%%% Figura ----------------------------------------------------------------
wid = 16;
hei = 7;
fig = figure('Units','centimeters','Position',[1 1 wid hei]);
N    = size(Y,2);
Clrs = colormap(winter(N));
LnSt = {'-','--','-.','-','--','-.'}; if N/length(LnSt)>1; LnSt = repmat(LnSt,1,ceil(N/length(LnSt))); end
%%% -----------------------------------------------------------------------
ax2 = axes;
pl1 = plot(f,Y);
for ii=1:N
    pl1(ii).Color     = Clrs(ii,:);
    pl1(ii).LineWidth = GraphProp.linewidth;
    pl1(ii).LineStyle = LnSt{ii};
end
set(ax2,GraphProp.Prop);set(ax2.XAxis,GraphProp.PropXA);set(ax2.YAxis,GraphProp.PropYA);set(ax2.Title,GraphProp.PropT);
set(ax2.XLabel,GraphProp.PropXL); set(ax2.YLabel,GraphProp.PropYL);
ax2.XLim = [0.1 MaxFreq];
ax2.Title.String  = '';
ax2.YLabel.String = YLabel;
ax2.XLabel.String = XLabel;
if nargin>5
    fig.Name = varargin{6};
end
%%% -----------------------------------------------------------------------
lgn1 = cellstr(num2str((1:N)'));
lgn1 = legend(lgn1,'FontSize',GraphProp.fontsize-2,'Location','northeast');
set(lgn1.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.9]));
end
%% PlotSpectrogram
function PlotSpectrogram(varargin)
%%%
y = varargin{1};
fs = varargin{2};
%%%
GraphProp = GraphicalProperties;
%%% Figura ----------------------------------------------------------------
wid = 16;
hei = 7;
for ii=1:size(y,2)
    fig = figure('Units','centimeters','Position',[1 1 wid hei],'Name',char(['Spectrogram: ',num2str(ii)]));
    spectrogram(y(:,ii),[],[],[],fs);
    ax1 = fig.Children(1);
    ax2 = fig.Children(2);
    set(ax1,'FontName',GraphProp.fontname,'FontSize',GraphProp.fontsize);
    set(ax2,GraphProp.Prop);
    set(ax2.XAxis,GraphProp.PropXA);
    set(ax2.YAxis,GraphProp.PropYA);
    set(ax2.Title,GraphProp.PropT);
    set(ax2.XLabel,GraphProp.PropXL);
    set(ax2.YLabel,GraphProp.PropYL);
    ax2.YLabel.String = 'Tiempo';
    ax2.XLabel.String = 'Frecuencia (Hz)';
    ax1.Label.String  = 'Potencia/frecuencia (dB/Hz)';
end
end
