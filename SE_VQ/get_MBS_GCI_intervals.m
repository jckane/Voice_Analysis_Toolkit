function interval = get_MBS_GCI_intervals(MBS,fs,T0mean,F0max)

% Function to detect intervals to be used for searching the LP-residual to
% determine GCI locations using the method described in Drugman et al
% (2012). 

%% Initial settings
if nargin < 4
    F0max=500;
end
F0max=F0max*2;
T0max=round(fs/F0max);
[~,idx]=findpeaks(MBS*-1,'minpeakdistance',T0max); % Find locations of negative peaks
N=length(idx);
search_rate=0.28;
search_left_rate=0.01;
interval=zeros(N,2);

%% Do processing
for n=1:N
   
    if length(T0mean)>1
        start=idx(n)-round(T0mean(idx(n))*search_left_rate);
        stop=idx(n)+round(T0mean(idx(n))*search_rate);
    else start=idx(n)-round(T0mean*search_left_rate);
        stop=idx(n)+round(T0mean*search_rate);
    end
    
    if start < 1
        start=1;
    end
    
    % Check start and end points of detected intervals
    if stop > length(MBS) && start < length(MBS)
        stop=length(MBS);
    elseif stop > length(MBS) && start >= length(MBS)
        break
    end
    
    interval(n,1)=start;
    interval(n,2)=stop;
end
    
    