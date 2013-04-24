function y = zeroPhaseHPFilt(x,fs,f_p,f_s,plots)

% Function to use a forward and backwards butterworth high pass filter to
% ensure zero phase 
%
% USAGE:    
%       Input:
%             x  : input signal
%             fs : sampling frequency (in Hz)
%             f_p : pass-band (in Hz)
%             f_s : stop band (in Hz)
%             plots : input 1 for a plot of the filters frequency response.
%
%       Output:
%             y  : filtered signal
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rp = 0.5;
Rs = 40;
Wp=f_p/(fs/2);
Ws=f_s/(fs/2);

[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[b,a]=butter(n,Wn,'high');

y = filtfilt(b,a,x);

if nargin > 4
    if plots==1
        freqz(b,a)
    end
end