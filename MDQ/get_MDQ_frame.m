function MDQ = get_MDQ_frame(res,fs,F0mean)

% Function to calculate the Maxima Dispersion Quotient (MDQ) from the
% inputted LP-residual of the speech signal. The method involves carrying 
% out zero-phase wavelet based filtering on the LP-residual of the input signal. 
% The dispersion of peaks across the different frequency bands are measured
% in relation to the peaks in the scale corresponding to the lowest frequency, 
% are averaged and then normalised to the local glottal period.
% Note that unlike the standard MDQ method, this particular approach does
% not require the extraction of GCIs, and may be more suitable for the
% analysis of speech recorded in less than ideal conditions. Informal
% testing has shown it to be broadly similar to the standard MDQ method but
% this has not been formally tested yet.
%
% INPUT:
%       res - Linear Prediction residual of speech signal
%       fs  - Sampling frequency (Hz)
%
% OUTPUT:
%       MDQ - Maxima Dispersion Quotient (values aligned to GCI)
%
%
% REFERENCE:
%       Kane, J., Gobl, C., ``Wavelet maxima dispersion for breathy to tense 
%       voice discrimination'', IEEE Trans. Audio Speech & Language
%       Processing [Under Review]


%% Initial settings
F0min=20;
F0max=500;
if nargin < 3
    [f0,VUV] = SRH_PitchTracking(x,fs,F0min,F0max);
    F0mean=median(f0(f0>F0min&f0<F0max&isnan(f0)==0&VUV==1));
end

T0mean=fs/F0mean;

i=2:6; % Analysis scales
s_num=length(i);

searchRate=0.2;
searchLen_cur=round(searchRate*T0mean); % Search region as a function of T0

MDQ=[];

winLen=round(32/1000*fs);
winShift=10/1000*fs;
start=1;
stop=start+winLen-1;

%% Do wavelet-based decomposition
[~,y_n] = do_aless_decomp(res,fs,i);

%% Do processing
cnt=1;
while stop <= length(res)
   
    low_y_cur=y_n(end-2,start:stop);
    [~,maxIdx]=max(low_y_cur);
    
    maxIdx_all=maxIdx+start-1;
    
    if maxIdx_all - searchLen_cur > 0
        start_frame=maxIdx_all - searchLen_cur;
        midpoint=searchLen_cur;
    else start_frame=1;
        midpoint=maxIdx_all;
    end
    if maxIdx_all + searchLen_cur <= length(res)
        stop_frame=maxIdx_all + searchLen_cur;
    else stop_frame=1;
    end
    
    y_cur=y_n(:,start_frame:stop_frame);
    
    dist_cur=zeros(1,s_num-1); % Allocate space for dispersion measure
    
    % Measure dispersion at each frequency band
    for m=1:s_num
        y_curn = y_cur(m,:);
        [~,maxIdx] = max(y_curn);
        dist_cur(m) = abs(midpoint-maxIdx);
    end  

    % Normalise average disperion to the local glottal period
   % MDQ(cnt)= std(dist_cur)/T0mean; 
    MDQ(cnt)= mean(dist_cur)/T0mean; 
    
    % Increment frames
    start=start+winShift-1;
    stop=start+winLen-1;
    cnt=cnt+1;
    
end
