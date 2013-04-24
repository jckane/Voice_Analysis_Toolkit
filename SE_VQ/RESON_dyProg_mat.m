function GCI_opt = RESON_dyProg_mat(GCI_relAmp,GCI_N,F0_mean,x,fs,trans_wgt,relAmp_wgt)

% Function to carry out dynamic programming method described in Ney (1989)
% and used previously in the ESPS GCI detection algorithm. The method
% considers target costs and transition costs which are accumulated in
% order to select the `cheapest' path, considering previous context

% USAGE: INPUT
%        GCI_relAmp - target cost matrix with N rows (GCI candidates) by M
%                     columns (mean based signal derived intervals).
%        GCI_N      - matrix containing N by M candidate GCI locations (in
%                     samples)
%        F0_inter   - F0 values updated every sample point       
%        x          - speech signal
%        fs         - sampling frequency
%        trans_wgt  - transition cost weight
%
%        OUTPUT
%        GCI        - estimated glottal closure instants (in samples)
% =========================================================================
% === FUNCTION CODED BY JOHN KANE AT THE PHONETICS LAB TRINITY COLLEGE ====
% === DUBLIN. 25TH October 2011 ===========================================
% =========================================================================

%% Initial settings
plots=0;
cost = (GCI_relAmp.*relAmp_wgt)';
ncands=size(GCI_N,1);
nframe=size(GCI_N,2);
GCI_N=GCI_N';
prev=zeros(nframe,ncands);      % traceback pointer
pulseLen = round(fs/F0_mean);
GCI_opt=zeros(1,nframe);

%% Do processing
for n=1:nframe
   
    if n>1
        costm=zeros(ncands,ncands);         % transition cost matrix: rows (previous), cols (current)
        
        for c=1:ncands
            % Transitions TO states in current frame
            start=GCI_N(n,c)-round(pulseLen/2);
            stop=GCI_N(n,c)+round(pulseLen/2);
            if stop > length(x)
                stop=length(x);
            end
            pulse_cur = x(start:stop);
            
            for p=1:ncands
                % Transitions FROM states in previous frame
                start=GCI_N(n-1,p)-round(pulseLen/2);
                stop=GCI_N(n-1,p)+round(pulseLen/2);
                if start<1
                    start=1;
                end
                if stop > length(x)
                    stop=length(x);
                end
                
                pulse_prev = x(start:stop);
                
                if isempty(pulse_cur) ||isnan(pulse_cur(1)) || ...
                        isnan(pulse_prev(1))
                    costm(p,c)=0;
                else
                    if length(pulse_cur)~=length(pulse_prev)
                        cor_cur=0;
                    else
                        cor_cur=corrcoef(pulse_cur(:),pulse_prev(:));
                        cor_cur=cor_cur(2);
                    end
                    costm(p,c) = (1-abs(cor_cur))*trans_wgt; % transition cost
                end
            end
        end
        
        costm=costm+repmat(cost(n-1,1:ncands)',1,ncands);  % add in cumulative costs
        [costi,previ]=min(costm,[],1);
        cost(n,1:ncands)=cost(n,1:ncands)+costi;
        prev(n,1:ncands)=previ;
    end
    
end

%% Do traceback
best=zeros(n,1);
[cbest,best(n)]=min(cost(n,1:ncands));
for i=n:-1:2
     best(i-1)=prev(i,best(i));
end

for n=1:nframe
    GCI_opt(n) = GCI_N(n,best(n));
end

%% Do plots
if plots
    GCI_norm=zeros(nframe,ncands);
    GCI_opt_norm=zeros(nframe,ncands);
    for n=1:nframe
        GCI_norm(n,:) = GCI_N(n,:)-GCI_N(n,1);
        GCI_opt_norm(n) = GCI_opt(n)-GCI_N(n,1);
    end
    subplot(211), plot(x), hold on, stem(GCI_N(:,1),ones(1,length(GCI_N(:,1)))*-.1,'r')
        stem(GCI_opt,ones(1,length(GCI_opt))*-.1,'k')
    subplot(212), plot(GCI_opt,GCI_norm,'rx'), ylim([-20 20]),
        hold on, plot(GCI_opt,GCI_opt_norm)
end