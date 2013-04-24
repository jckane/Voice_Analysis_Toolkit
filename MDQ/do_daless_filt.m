function y_i = do_daless_filt(s,h_i)

% Function to convolve a given input signal, s, with the wavelet function
% h_i

s_h_i=conv(s,h_i);
halfLen = ceil(length(h_i)/2);
y_i = s_h_i(halfLen:halfLen+length(s)-1);


