function GCI = align_GCIs_to_residual_peaks(GCI_in,res,num_peaks,search_rate)

%% Initial settings
if nargin < 3
    num_peaks=5;
end
if nargin < 4
    search_rate=.2;
end
N=length(GCI_in);
GCI=zeros(N,num_peaks);
res_amp=zeros(N,num_peaks);
GCI(end,:)=GCI_in(end,:);

%% Do processing
for n=1:N-1
    start=GCI_in(n);
    peak_next=GCI_in(n+1);
    T0=peak_next-start+1;
    
    stop=round(start+(T0*search_rate));
    res_frame=res(start:stop);
    [amp,idx]=sort(res_frame,'descend');
    res_amp(n,:)=amp(1:num_peaks);
    GCI(n,:)=idx(1:num_peaks)+start-1;
    
end
