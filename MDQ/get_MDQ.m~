function MDQ = get_MDQ(res,fs,GCI)

% Function to calculate the Maxima Dispersion Quotient (MDQ) from the
% inputted LP-residual of the speech signal. The method involves carrying 
% out zero-phase wavelet based filtering on the LP-residual of the input signal. 
% The dispersion of peaks across the different frequency bands are measured
% in relation to the glottal closure instant, are averaged and then
% normalised to the local glottal period
%
% INPUT:
%       res - Linear Prediction residual of speech signal
%       fs  - Sampling frequency (Hz)
%       GCI - Glottal closure instants (in samples)
%
% OUTPUT:
%       MDQ - Maxima Dispersion Quotient (values aligned to GCI)
%
%
% REFERENCE:
%       Kane, J., Gobl, C., (2013)``Wavelet maxima dispersion for breathy to tense 
%       voice discrimination'', IEEE Trans. Audio Speech & Language
%       Processing, 21(6), pp. 1170-1179.
%
% NOTES:
%       Performance appears to be optimal for F0 in the range [50, 200],
%       between [200, 300] it is satisfactory but after 300 Hz the
%       performance deteriorates significantly. One idea would be to
%       measure the strength of periodicity in the higher scales (lower
%       frequencies) to determine whether higher ones should be omitted
%       from the dispersion measurement


%% Initial settings
GCI=unique(GCI);

i=2:6; % Analysis scales [2000, 1000, 500, 250, 125]
s_num=length(i);

searchRate=0.2;

MDQ=zeros(1,length(GCI));

%% Do wavelet-based decomposition
[~,y_n] = do_aless_decomp(res,fs,i);

%% Do processing
for n=1:length(GCI)
   
    % Get local Glottal period
    if n==1
        T0=GCI(n+1)-GCI(n);
    else T0=GCI(n)-GCI(n-1);
    end
    
    searchLen_cur=round(searchRate*T0); % Search region as a function of T0
    
    
    if GCI(n)-T0 > 0 && GCI(n)+T0 <= length(res)

        % Get refined frame start and end points
        start_ser=GCI(n)-searchLen_cur;
        finish_ser=GCI(n)+searchLen_cur;
        midpoint=searchLen_cur;
        
        % Current frame
        y_cur=y_n(:,start_ser:finish_ser);

        dist_cur=zeros(1,s_num); % Allocate space for dispersion measure

        % Measure dispersion at each frequency band
        for m=1:s_num
            y_curn = y_cur(m,:);
            [~,maxIdx] = max(y_curn);
            dist_cur(m) = abs(midpoint-maxIdx);
        end  

        % Normalise average disperion to the local glottal period
        MDQ(n)= mean(dist_cur)/T0; 
       
     end
end

MDQ(isinf(MDQ))=0; % Remove any infinity values

