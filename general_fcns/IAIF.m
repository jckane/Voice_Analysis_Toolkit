function [g_iaif,ar_lpc,e_lpc] = IAIF(x,fs,GCI,p)

% Function to carry out iterative and adaptive inverse filtering (Alku et
% al 1992). 

% USAGE:    
%       Input:
%             x   : speech signal (in samples)
%             fs  : sampling frequency (in Hz)
%             GCI : Glottal closure instants (in samples)
%             p   : LPC prediction order
%
%       Output:
%             g_iaif : glottal flow derivative estimate
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initial settings
x=x(:);

if nargin < 4
    p = round(fs/1000)+2; % LPC filter order
end

if nargin < 3,
    GCI = SE_VQ(x,fs); 
end


% ------------------------------------------------
% High-pass filter to eliminate DC component - Part 1
nc=704; %352
fc=50; % 30 Hz in Alku (1992)
nfc=fc/(fs/2);
b_hp=fir1(nc,nfc,'high');
[x_hp]=filter(b_hp,1,x);
x=x(:)';
x_hp=x_hp(:)';
x_filt=[x_hp(nc/2+1:end) x(end-nc/2+1:end)];

% ------------------------------------------------
% emphasise high-frequencies of speech signal (LPC order 1) - PART 2 & 3
ord_lpc1=1;
[x_emph,ar_lpc1,e_lpc1]=calc_residual(x_filt,x_filt,ord_lpc1,GCI);

% ------------------------------------------------
% first estimation of the glottal source derivative - PART 4 & 5
ord_lpc2=p;
[residual1,ar_lpc2,e_lpc2]=calc_residual(x_filt,x_emph,ord_lpc2,GCI);

% integration of the glottal source derivative to calculate the glottal
% source pulse - PART 6 (cancelling lip radiation)
ug1=residual1;
% ------------------------------------------------
% elimination of the source effect from the speech spectrum - PART 7 & 8

ord_lpc3=4;
[vt_signal,ar_lpc3,e_lpc3]=calc_residual(x_filt,ug1,ord_lpc3,GCI);

% ------------------------------------------------
% second estimation of the glottal source signal - PART 9 & 10

ord_lpc4=p;
[residual2,ar_lpc,e_lpc]=calc_residual(x_filt,vt_signal,ord_lpc4,GCI);

% -----------------------------------------------
% upsampling of the residual to the original sampling rate

g_iaif=residual2;


