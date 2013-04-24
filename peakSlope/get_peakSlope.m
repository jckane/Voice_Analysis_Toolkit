function peakSlope = get_peakSlope(s,fs)

% Function to calculate peakSlope measurement on a fixed frame basis.
% Unlike in the Interspeech paper maxima are converted to dB.

% INPUT
%        s    - speech signal in samples
%        fs   - sampling frequency       
%
% OUTPUT
%        peakSlope - peakSlope parameter (measured every 10 ms)
%
% REFERENCE
%        Kane, J., Gobl, C., (2011) ``Identifying regions of non-modal phonation 
%        using features of the wavelet transform'', Proceedings of
%        Interspeech
%
% NOTES
%        The dependency on signal energy has now been address by modifying
%        the measurement described in the publication and now fitting the
%        regression line to log10 of the maxima.
%
% =========================================================================
% === FUNCTION CODED BY JOHN KANE AT THE PHONETICS LAB TRINITY COLLEGE ====
% === DUBLIN. 8th June 2011================================================
% =========================================================================

%% Initial settings
frameLen_ms = 40; % Frame length chosen to ensure one pulse length down to f0=25 Hz
frameShift_ms = 10; % Frame shift set to 10 ms
frameLen = (frameLen_ms/1000)*fs; % Convert frame length to samples
frameShift = (frameShift_ms/1000)*fs; % Convert frame shift to samples

peakSlope=zeros(1,round((length(s)-frameLen)/frameShift));

%% Do wavelet decomposition
i=0:6; % i=0:6 => 8 kHz, 4 kHz, 2 kHz, 1 kHz, 500 Hz, 250 Hz, 125 Hz

y=zeros(length(i),length(s)); % Allocate space for the different frequency bands

for n=1:length(i)
    
    h_i = daless_MW(i,n,fs); % Generate mother wavelet
    y_i = do_daless_filt(s,h_i); % Carry out zero-phase filtering
    y(n,:) = y_i;
end

%% Measure peakSlope per frame
start=1;
finish = start+frameLen-1;
m=1;

while finish <= length(s)
    maxima = zeros(1,length(i)); % allocate space
    
    for n=1:length(i)
        maxima(n) = max(abs(y(n,start:finish))); % measure peaks at each scale
    end
    maxima = log10(maxima(end:-1:1)); % reverse order to follow frequency order and convert to dB
    t=1:length(maxima);
    p=polyfit(t,maxima,1); % do straight line regression fitting
    peakSlope(m) = p(1); % take slope coefficient from regression line
    
    m=m+1;
    start = start+frameShift;
    finish = start+frameLen-1;
end
peakSlope(isnan(peakSlope)==1)=0;

