function y = zeroPhaseLPFilt(x,fs,f_p,f_s,plots)

% Function to use a forward and backwards butterworth high pass filter to
% ensure zero phase 

Rp = 0.5;
Rs = 6;
Wp=f_p/(fs/2);
Ws=f_s/(fs/2);

[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[b,a]=butter(n,Wn,'low');

y = filtfilt(b,a,x);

if nargin > 4
    if plots==1
        freqz(b,a)
    end
end