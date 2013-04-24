function [y,y_norm] = do_aless_decomp(s,fs,i)

% Function to do full wavelet based decomposition of an input signal

if nargin < 3
    i=0:6; % i=0:6 => 8 kHz, 4 kHz, 2 kHz, 1 kHz, 500 Hz, 250 Hz, 125 Hz
end

%% Allocate space
y=zeros(length(i),length(s));
y_norm=zeros(length(i),length(s));

%% Do filtering
for n=1:length(i)
    
    h_i = daless_MW(i,n,fs);
    y_i = do_daless_filt(s,h_i);
    y(n,:) = y_i;
    y_norm(n,:) = y_i/max(abs(y_i));
end