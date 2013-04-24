function [H2H1,res_p,ZCR,F0,F0mean,enerN,pow_std,creakF0] = get_KD_creak_features(x,fs)

% Function to extract the various features used for creak detection by the
% binary decision tree classifier
%
% USAGE:    
%       Input:
%             x   : speech signal
%             fs  : sampling frequency (in Hz)
%             res : Linear Prediction residual of speech signal
%
%       Output:
%             creak_dec: binary creak decision vector
%
% REFERENCE:
%       Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%
%       Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', 27(4), pp. 1028-1047, Computer Speech and Language.
%
%       Drugman, T., Kane, J., Gobl, C., (2012) `Resonator-based creaky 
%       voice detection', Interspeech 2012, Portland, Oregon, USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial settings
F0min=20;
F0max=500;
winShift=round(10/1000*fs); % Sample every 10 ms
res = GetLPCresidual(x,round(25/1000*fs),round(5/1000*fs),round(fs/1000)+2);
[f0,VUV,~] = SRH_PitchTracking(x,fs,F0min,F0max);
F0mean=median(f0(VUV==1&f0>F0min&f0<F0max));

%% Extract features for creak detection
time=winShift:winShift:length(x);

[~,Es,ZC_ms,Xpos] = sil_unv_features(x,fs,32);
ZC_inter = interp1(Xpos,ZC_ms,1:length(x));
ener_norm = Es-max(Es);
Es_inter=interp1(linspace(1,length(x),length(ener_norm)),ener_norm,1:length(x));
[~,~,pow_std_inter] = get_short_pow(x,fs);


[H2H1,F0_creak] = get_creak_H2H1(res,fs,F0mean);
peak_inter = res_peak(x,fs,F0mean,res,Es);

%% Do resampling
Es_re=Es_inter(time);
pow_re=pow_std_inter(time);
H2H1_re=H2H1(time);
peak_re=peak_inter(time);
peak_re(isnan(peak_re))=0;
ZC_re=ZC_inter(time);
ZC_re(isnan(ZC_re))=0;
F0_creak=F0_creak(time);

%% Save to structs
ZCR=ZC_re(:);
enerN=Es_re(:);
creakF0=F0_creak(:);
pow_std=pow_re(:);
H2H1=H2H1_re(:);
res_p=peak_re(:);
F0mean=zeros(size(res_p))+F0mean;
F0=interp1(linspace(1,length(F0mean),length(f0)),f0,1:length(F0mean));




