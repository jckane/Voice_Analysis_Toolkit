function [GCI,rep,res,MBS] = SE_VQ(x,fs,f0,VUV,creak)

% Function to extract GCIs using an adapted version of the SEDREAMS
% algorithm which is optimised for non-modal voice qualities. Ncand maximum
% peaks are selected from the LP-residual signal in the interval defined by
% the mean-based signal. A dynamic programming algorithm is then used to
% select the optimal path of GCI locations. Then a post-processing method,
% using the output of a resonator applied to the residual signal, is
% carried out to remove false positives occurring in creaky speech regions
%
% Usage:
% 
% Input: 
%           x            - speech signal (in samples)
%           fs           - sampling frequency (Hz)
%           f0           - fundamental frequency contour (Hz) - OPTIONAL
%           VUVDecisions - Voiced/unvoiced decision (binary) - OPTIONAL
%           creak        - creaky voice binary decision - OPTIONAL
%
%
% REFERENCE:
%       Kane, J., Gobl, C., (2013) `Evaluation of glottal closure instant 
%       detection in a range of voice qualities', Speech Communication 
%       55(2), pp. 295-314.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION WAS CODED BY JOHN KANE AT THE PHONETICS AND SPEECH LAB IN %%%%%
%% TRINITY COLLEGE DUBLIN ON 2013. THE SEDREAMS FUNCTION WAS%%%%%
%% CODED BY THOMAS DRUGMAN OF THE UNIVERSITY OF MONS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
F0min=20;
F0max=500;

if nargin < 3 || isempty(f0) || isempty(VUV)
    % f0 and voicing decision is not critical and this can be replaced with
    % a different algorithm is necessary
    [f0,VUV] = SRH_PitchTracking(x,fs,F0min,F0max);  
end
F0mean=median(f0(f0>F0min&f0<F0max&VUV==1));
F0max=max(medfilt1(f0(VUV==1),13));
if F0mean < 70
    disp('Utterance likely to contain creak')
    F0mean=80;
end
T0mean = fs/F0mean; % Rough period length for mean-based signal   

winLen = 25; % window length in ms 
winShift = 5; % window shift in ms 
LPC_ord = (fs/1000)+2; % LPC order
Ncand=5; % Number of candidate GCI residual peaks to be considered in the dynamic programming

trans_wgt=1; % Transition cost weight
relAmp_wgt=0.3; % Local cost weight

repNum=2;
removeThresh=0.4; % Threshold for removing false GCIs
search_reg=1.3/1000*fs;

%% Calculate LP-residual and extract N maxima per mean-based signal determined intervals
res = GetLPCresidual(x,winLen/1000*fs,winShift/1000*fs,LPC_ord); % Get LP residual
rep = RCVD_reson_GCI(res,fs,F0mean); % Get resonator output
MBS = get_MBS(x,fs,T0mean); % Extract mean based signal
interval = get_MBS_GCI_intervals(MBS,fs,T0mean,F0max); % Define search intervals
[GCI_N,GCI_relAmp] = search_res_interval_peaks(res,interval,Ncand); % Find residual peaks
GCI = RESON_dyProg_mat(GCI_relAmp',GCI_N',F0mean,x,fs,trans_wgt,relAmp_wgt); % Do dynamic programming

%% Remove false alarms as weak peaks in resonator output
if nargin > 4 && length(creak)==length(x)
    disp('Doing post-processing in detected creaky voice regions')
    GCI = GCI_creak_postproc(GCI,creak,search_reg,rep,removeThresh,repNum);
end

