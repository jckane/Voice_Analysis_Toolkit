function [GCI,rel_amp] = search_res_interval_peaks(res,interval,Ncand)

% Function to search for Ncand peaks in the residual within each search
% interval
N=size(interval,1);
GCI=zeros(N,Ncand);
rel_amp=zeros(N,Ncand);

for n=1:N
   
    start=interval(n,1);
    stop=interval(n,2);
    
    if stop <= start
        GCI(n,:)=NaN;
    elseif stop-start < Ncand
        [amp,idx]=max(res(start:stop));
        GCI_cur=GCI_cur+start-1;
        
        GCI(n,:)=GCI_cur;
        rel_amp(n,:)=0;
        
    else [amp,idx]=sort(res(start:stop),'descend');
        GCI_cur=idx(1:Ncand)+start-1;
        GCI(n,:)=GCI_cur(:)';
        rel_amp(n,:)=1-(amp(1:Ncand)/max(amp(1:Ncand)));
    end
end

