function h_i = daless_MW(i,n,fs)

% Function to generate the mother wavelet from d'Alessandro 2011.

s=2.^i; % convert to octave bands
f_o = fs/2;
tau = 1/(2*f_o);

t=(-1000:1000)./fs;
mother_wavelet = @(s_n) -cos(2*pi*f_o.*(t./s_n)).*exp(-((t./s_n).^2)/(2*(tau^2))); % mother wavelet
h_i = mother_wavelet(s(n));