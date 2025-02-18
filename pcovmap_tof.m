% covariance mapping for a the 20pC Ne tof measurement in LCLS
% based on L. Frasinski work from 2011
close all
clear all
clc
% Partial covariance expression:
% pcovYX = covYX - covIY * covIX / varI

% where :
% covYX = <(Y - avY)*(X - avX)> = avYX - avY * avX is the covariance of the spectra (matrix)
% covIX = <(I - avI)*(X - avX)> = avIX - avI * avX is the covariance of the intensity spectrum (row vector)
% covIY = <(I - avI)*(Y - avY)> = avIY - avI * avY is the covariance of the intensity spectrum (col vector)
% varI  = <(I - avI)^2> = avI2 - avI^2 is the variance of intensity (scalar)
% avYX  = <Y*X> is the avg spectrum-spectrum product (matrix)
% avIX  = <I*X> is the avg intensity-spectrum product (row vector)
% avIY  = <I*Y> is the avg intensity-spectrum product (col vector)
% avI2  = <I^2> is the avg intensity squared (scalar)
% avX   = <X>   is the avg spectrum (row vector)
% avY   = <Y>   is the avg spectrum (col vector)
% avI   = <I>   is the avg intensity (scalar)

% and,
% <...> is averaging over n
% n = 1 : nMax           laser shots
% i = 1 : iMax           tof samples
% j = 1 : jMax           energy samples
% X = X(n,i) or X(n,j)   rows of spectra
% Y = X(i)' or X(j)'     cols of spectra
% I = I(n)               shot intensities, col vector


% % Analysis run parameters
%runName1 = 'e88-r0105';
%runName2 = '';               % = '' if not combining data sets
%fMax1 = 103;                 % number of files in the run
%fMax2 =   0;                 % = 0 if not combining data sets
%sample = 'Ne, 20pC, 0V';     % used in figure titles
%U_retarding = 0;             % retarding potential /V
%E = linspace(U_retarding,250, 600); % energy scale /eV
ff = 1;                       % fudge factor for intensity correction

%fMax = fMax1 + fMax2;
% % Magnetic bottle parameters by Timur Osipov
% %     Uret/V E0/eV T0/us L0/cm
% ETL = [  0   0.110 3.926 218.303
%        100 102.704 3.931 211.604
%        150 155.345 3.937 205.866];
% ind = find(ETL(:,1) == U_retarding, 1);
% if isempty(ind)
%   error(['Missing T0 and L0 for retarding potential ' num2str(U_retarding) ' V']);
% end
% E0 = ETL(ind,2);            % energy reduction due to U_retarding
% T0 = ETL(ind,3);            % zero tof us
% L0 = ETL(ind,4);            % length of tof tube/cm

%% Load or obtain covariance map parameters
if exist('CovMap105.mat', 'file')
    load('CovMap105.mat');

else % what's under the hood to get covariance parameters from data (starting row 121):
    covName=['CovMap' num2str(runName1(end-2:end))];
    % read iMax and eMax
    dataName = [runName1 '-001.mat'];
    dat = load(dataName, 'eMax','iMax');
    iMax = dat.iMax;
    nMax0 = fMax*dat.eMax; 

    % read data
    textprogressbar('get cov parameters: '); tic

    n = 0;
    for fn = 1:fMax
        % combine runs
        if fn <= fMax1
            fStr = sprintf('%03d', fn);     % string of file number
            dataName = [runName1 '-' fStr '.mat'];
        else
            fStr = sprintf('%03d', fn-fMax1);
            dataName = [runName2 '-' fStr '.mat'];
        end

        % read common data and intensity
        if fn == 1
            dat = load(dataName);
            W0 = mean(mean(dat.W(:, end-1000:end))); % zero level assumption 
            T      = 1E6*dat.T;      % in us
            gain   = 1E3*dat.gain;   % in mV
            offset = 1E3*dat.offset; % in mV
        else
            dat = load(dataName, 'eMax','GasDetector','W');
        end
        dat.I = mean(dat.GasDetector(:,3:4), 2); % after attenuation

        % remove events with zero intensity
        ee = (dat.I == 0);        % zero-intensity indices
        dat.I(ee)   = [];         % remove intensities
        dat.W(ee,:) = [];         % remove waveforms
        dat.eMax = length(dat.I); % count only good events

        % tof to energy
        [dat.X, E] = t2e(T,gain*(double(dat.W) - W0),E, L0,T0,E0);  % in mV*us/eV,
        jMax = length(E);

        % we kept some of single-shot spectra (later plotted)
        W = dat.W(1:20, :);
        X = dat.X(1:20, :);
        
        if fn<2
            % initialize arrays
            avYX = zeros(jMax,jMax);
            I    = zeros(nMax0,1);
            Wsum = zeros(nMax0,1);
            avU  = zeros(1,iMax);
            avX  = zeros(1,jMax);
            avIX = zeros(1,jMax);
            avI2 = 0;
            avI  = 0;
        end

        % *** get covariance parameters ***
        for en = 1 : dat.eMax
            x = dat.X(en,:);
            avYX = avYX + x' * x;
            avIX = avIX + dat.I(en) * x;
        end

        avX  = avX  + sum(dat.X);
        avU  = avU  + gain*sum(double(dat.W) - W0);
        avI  = avI  + sum(dat.I);
        avI2 = avI2 + sum(dat.I .^ 2);
        I(n+1 : n+dat.eMax) = dat.I;
        Wsum(n+1 : n+dat.eMax) = sum(double(dat.W)-W0, 2);
        n = n + dat.eMax;

        textprogressbar((fn-1)/fMax*100);
    end

    nMax = n;
    avYX = avYX/nMax;
    avIX = avIX/nMax;
    avU  = avU/nMax;
    avX  = avX/nMax;
    avI  = avI/nMax;
    avI2 = avI2/nMax;
    nn = (I==0);   % ind of empty shots
    I(nn)    = [];
    Wsum(nn) = [];
    clear n nn en x fn fStr dataName dat
    save(covName,'avYX','avIX','avU','avX','avI','avI2','I','Wsum','T','E','X','W','W0', 'gain','offset','nMax0','nMax','iMax','jMax','runName1');

    textprogressbar(['done! ' '(' num2str(toc) ' sec)']);
end

%% get partial covariance
avY = avX';
covIX = avIX - avI * avX;
varI  = avI2 - avI^2;
covYX  = (avYX - avY * avX);        % covariance
Icor = covIX' * covIX / varI;       % decorrelation of intensity
pcovYX = covYX - ff * Icor;         % partial covariance

%% plot it all...
figure('Position',[50 50 700 700])

% inspecting some single shots...
subplot(2,2,1)
cmap=jet(size(W,1));
for nP=1:size(W,1)
    plot(T, 20*(nP-1)+gain*(double(W(nP,:)) - W0),'Color',cmap(nP,:));
    hold on
end
xlabel('tof (\mus)')
ylabel('adc signal (mV)')
xlim([4 5])
ylim([-10 max(20*(nP-1)+gain*(double(W(nP,:)) - W0))])
title('Single-shot spectra')

subplot(2,2,2)
for nP=1:size(W,1)
    plot(E, 20*(nP-1)+X(nP,:),'Color',cmap(nP,:));
    hold on
end
xlabel('energy (eV)')
ylabel('converted signal (mV\mus/eV)')
ylim([-10 max(20*(nP-1)+X(nP,:))])

% avg tof spectrum
subplot(2,2,3)
plot(T, avU)
xlabel('tof  (\mus)')
ylabel('adc signal (mV)')
title('avg spectra');
xlim([4 5])

% in energy...
subplot(2,2,4)
plot(E, avX)
xlabel('energy (eV)')
ylabel('converted signal  (mV\mus/eV)')

hcb = colorbar;
hcb.Position = [0.94 0.11 0.02 0.815];  % Adjust these values as needed
colormap(jet(size(W,1)))
caxis([1 size(W,1)+1])
hcb.Limits = [1, size(W,1)+1];
hcb.Ticks = 0.5+(1:size(W,1));
hcb.TickLabels = arrayfun(@num2str, 1:size(W,1), 'UniformOutput', false);
title(hcb, 'shot #');

%% Intensity correction covariance
figure('Position',[200 200 300 300])
imagesc(E, E, Icor)
set(gca,'Ydir','normal')
xlabel('energy  (eV)')
ylabel('energy  (eV)')
title('Intensity correction (Icor)')
caxis([0 1])

%% partial covariance map
figure
% axis limits and ranges
xL = [min(E) max(E)];
xR = xL(2) - xL(1);
yL = [0 1.1*max(avX)];
yR = yL(2) - yL(1);

% Partial covariance
subplot(3, 3, [2, 3, 5, 6]);  % Combine four subplot locations into one larger plot
imagesc(E, E, pcovYX)
set(gca,'Ydir','normal')
caxis([0 2])
title('Partial covariance  (mV\mus/eV)^2')


% bottom tof spectra
subplot(3, 3, [8, 9]);  % Last row, last three columns
plot(E,avX,'r')
ylabel('signal (mV\mus/eV)')
xlabel('energy  (eV)')
xlim(xL)
ylim(yL)

% left tof spectrum
subplot(3, 3, [1, 4]);  % Span 3 rows in the 1st column of a 3x4 grid
plot(E,avX,'r')
ylabel('signal (mV\mus/eV)')
xlabel('energy  (eV)')
xlim(xL)
ylim(yL)
view(-90,90)

% Aux functions
function [Y,E] = t2e(T,X,E, L0,T0,E0)
% Converts tof to energy
% Inputs:
%   T(i)   - tof scale (us), row vector
%   X(n,i) - rows of tof spectra
%   E(j)   - energy scale (eV)
% Conversion parameters:
%   L0     - length of tof tube (cm)
%   T0     - optional zero tof (us) = 0 if omitted
%   E0     - optional retarding potential (V) = 0 if omitted
% Outputs:
%   Y(n,k) - rows of energy spectra (Xunits*ns/eV)
%   E(k)   - energy scale /eV, same as E(j) but any E outside the range of T is rejected

%  Default parameters
if nargin < 6
    E0 = 0;
end
if nargin < 5
    T0 = 0;
end

% Check inputs
if ~isvector(T)
    error('T must be a vector')
end
if T(1) > T(end)
    error('T must be increasing')
end
iMax = length(T);      % number of tof samples
if ~isvector(E)
    error('E must be a vector')
end
if E(1) > E(end)
    error('E must be increasing')
end
jMax = length(E);      % initial number of energy points
if iscolumn(T)
    T = T';
end
if iscolumn(E)
    E = E';
end
if size(X,2) ~= iMax
    error('Length of T must match number of cols in X')
end
nMax = size(X,1);      % number of spectra

% Find integration limits
E1 = ([2*E(1)-E(2) E] + [E 2*E(end)-E(end-1)])/2; % energy mid points
E1(E1<E0) = [];                  % remove negative energies
T1 = tof2(E1, L0,T0,E0);         % tof integration starting points
jj = T1 < T(1) | T1 > T(end);    % indices of T1 outside T range
T1(jj) = [];                     % remove T1 outside T range
E1(jj) = [];                     % remove corresponding E1
X1 = interp1(T',X', T1')';       % find tof signal interpolated at T1
E = (E1(1:end-1) + E1(2:end))/2; % new energy scale
kMax = length(E);                % to label energy use index k from now

% Convert energy point by energy point
Y = zeros(nMax, kMax);               % preallocate output array
for k = 1 : kMax
    ii = T1(k+1) < T & T < T1(k);      % indices between integration limits
    T12 = [T1(k+1)   T(ii)   T1(k)];   % integration tofs
    X12 = [X1(:,k+1) X(:,ii) X1(:,k)]; % spectrum parts to integrate
    Y(:,k) = trapz(1E3*T12, X12, 2);   % integrate in ns and store
end

end

function textprogressbar(c)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m
% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version
% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/
%% Initialization
persistent strCR;           %   Carriage return pesistent variable
% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 20;   %   The total number of dots in a progress bar
%% Main
if isempty(strCR) && ~ischar(c),
    % Progress bar must be initialized with a string
    %error('The text progress must be initialized with a string');
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;

elseif isempty(strCR) && ischar(c),
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c),
    % Progress bar  - termination
    strCR = [];
    fprintf([c '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];

    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end

    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);

else
    % Any other unexpected input
    error('Unsupported argument type');
end
end