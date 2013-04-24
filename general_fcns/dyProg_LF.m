function [params,cost,best,Rd_n] = dyProg_LF(glot,fs,GCI,plots)

% Function to fit LF models to glottal source pulses by using an exhaustive
% search method and subsequent dynamic programming, using a similar
% approach to that in RAPT (Talkin 1990). A final optimisation procedure is
% carried out to refine the fit by varying Ra Rk and Rg
%
% For further details see: 
%       Kane, J., Gobl, C. (2013) ``Automating manual user strategies for
%       precise voice source analysis'', Speech Communication 55(3), pp.
%       397-414.
%
%       AND
%
%       Kane, J., Yanushevskaya, I., Ní Chasaide, A.,
%       Gobl, C., (2012) ``Exploiting time and frequency domain measures for precise
%       voice source parameterisation'', Proceedings of Speech Prosody, 
%       Shanghai, China
%
% Usage:
%       Input:
%             glot  - glottal source (derivative) waveform in samples
%             fs    - sampling frequency (Hz)
%             GCI   - glottal closure instants (samples)
%             plots - 0 or 1 to indicate whether to output plot (optional)
%
%       Output:
%             params - struct containing LF model based parameters {F0,EE,Ra,Rk,Rg,OQ}
%             cost   - cummulated cost matrix from dynamic programming
%             best   - indices of optimal Rd values
%             Rd_n   - matrix of Rd candidates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coded by John Kane in the Phonetics and Speech Lab, Trinity College,
%% Dublin, June 2012
%% Updated February 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial settings
MVF=2500; % Maximum voiced frequency set to 2500 Hz
glot=glot(:);

% Dynamic programming weights
time_wgt=0.1;
freq_wgt=0.3;
trans_wgt=0.3;

% Allocated space
params.F0=zeros(1,length(GCI));
params.EE=zeros(1,length(GCI));
params.Ra=zeros(1,length(GCI));
params.Rk=zeros(1,length(GCI));
params.Rg=zeros(1,length(GCI));
params.UP=zeros(1,length(GCI));

options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',20);
EE=zeros(1,length(GCI));
Rd_set=[0.3:0.17:2];
pulseNum=2;

% Dynamic programming settings
nframe=length(GCI);
ncands = 5; % Number of candidate LF model configurations to consider
Rd_n=zeros(nframe,ncands);  
cost=zeros(nframe,ncands);      % cumulative cost
prev=zeros(nframe,ncands);      % traceback pointer

GCI = correct_GCI_glot(glot,GCI,fs); % align GCIs to main excitations

ss = get_spec_stat(glot,fs,GCI); % Get spectral similiarity values
ss = smooth(medfilt1(ss,5),5);
% flag=0;

%% Do processing - exhaustive search and dynamic programming
for n=1:length(GCI)
    
    % get framing information
    if n==1
        pulseLen=round((GCI(n+1)-GCI(n))*pulseNum);
        F0_cur=fs/(round(GCI(n+1)-GCI(n)));
    else pulseLen=round((GCI(n)-GCI(n-1))*pulseNum);
        F0_cur=fs/(round(GCI(n)-GCI(n-1)));
    end
    pulseLen=abs(pulseLen);
    
     if GCI(n)-round(pulseLen/2) > 0
        start=GCI(n)-round(pulseLen/2);
    else start=1;
    end
    
    if GCI(n)+round(pulseLen/2) <= length(glot)
        finish = GCI(n)+round(pulseLen/2);
    else finish = length(glot);
    end

  %  [start finish]
    % Get fft spectrum
%     if finish-start > 0 
        glot_seg=glot(start:finish).*hanning(finish-start+1);
        glot_seg=glot_seg(:);
        glot_seg_spec=20*log10(abs(fft(glot_seg)));
        freq=linspace(0,fs,length(glot_seg));

        err_mat=zeros(1,length(Rd_set));
        err_mat_time=zeros(1,length(Rd_set));
        EE(n)=abs(min(glot_seg));


        %% exhaustive search
        for m=1:length(Rd_set)
            [Ra_cur,Rk_cur,Rg_cur] = Rd2R(Rd_set(m),EE(n),F0_cur);

            % Generate test LF pulse and get spectral information
            pulse = lf_cont(F0_cur,fs,Ra_cur,Rk_cur,Rg_cur,EE(n));
            if isnan(pulse(1))
                err_mat(1:length(Rd_set))=1;
                break
            end

            LFgroup = makePulseCentGCI(pulse,pulseLen,GCI(n)-start,finish-GCI(n));
            LFgroup_win=LFgroup(:).*hanning(finish-start+1);
            LFgroup_win_spec=20*log10(abs(fft(LFgroup_win)));

            LFgroup_win=LFgroup_win(:);

            % Time domain error function
            cor_time = corrcoef(glot_seg,LFgroup_win);
            cor_time=abs(cor_time(2));
            err_time=1-cor_time;
            err_mat_time(m)=err_time;

            % Frequency domain error function
            cor_freq = corrcoef(glot_seg_spec(freq<MVF),LFgroup_win_spec(freq<MVF));
            cor_freq=abs(cor_freq(2));
            err_freq=1-cor_freq;

            % Combined error with weights
            err_mat(m)=(err_time*time_wgt)+(err_freq*freq_wgt);

        end

        % Find best ncands (local costs and Rd values)
        [err_mat_sort,err_mat_sortIdx]=sort(err_mat);
        Rd_n(n,1:ncands)=Rd_set(err_mat_sortIdx(1:ncands));
        exh_err_n=err_mat_sort(1:ncands); 
        cost(n,1:ncands) = exh_err_n(:)';

        %% Find optimum Rd value (dynamic programming)
        if n>1
            costm=zeros(ncands,ncands);         % transition cost matrix: rows (previous), cols (current)

            for c=1:ncands
                % Transitions TO states in current frame
                [Ra_try,Rk_try,Rg_try] = Rd2R(Rd_n(n,c),EE(n),F0_cur);
                LFpulse_cur = lf_cont(F0_cur,fs,Ra_try,Rk_try,Rg_try,EE(n));

                for p=1:ncands
                    % Transitions FROM states in previous frame
%                     if flag==1
%                         Rd_n(n-1,p)=median(Rd_n(1:n-2,p));
%                     end
                    [Ra_prev,Rk_prev,Rg_prev] = Rd2R(Rd_n(n-1,p),EE(n),F0_cur);
                    LFpulse_prev = lf_cont(F0_cur,fs,Ra_prev,Rk_prev,Rg_prev,EE(n));

                    if isnan(LFpulse_cur(1)) || isnan(LFpulse_prev(1))
                        costm(p,c)=0;
                    else
                        cor_cur=corrcoef(LFpulse_cur(:),LFpulse_prev(:));
                        cor_cur=cor_cur(2);
                        costm(p,c) = (1-abs(cor_cur))*trans_wgt; % transition cost
                    end
                end
                costm=costm.*ss(n); % Use spectral similarity measure to reduce transition effect in when there are rapid changes
            end
            costm=costm+repmat(cost(n-1,1:ncands)',1,ncands);  % add in cumulative costs
            [costi,previ]=min(costm,[],1);
            cost(n,1:ncands)=cost(n,1:ncands)+costi;
            prev(n,1:ncands)=previ;
        end
%         flag=0;
%     else flag=1;
%     end
end
    
%% Do traceback
best=zeros(n,1);
[~,best(n)]=min(cost(n,1:ncands));

for i=n:-1:2
     best(i-1)=prev(i,best(i));
end

Rd_opt=zeros(1,nframe);
for n=1:nframe
%     if n > length(best)
%         break
%     end
    Rd_opt(n) = Rd_n(n,best(n));
end

Rd_opt = smooth(medfilt1(Rd_opt,11),5)*.5;
params.Rd=Rd_opt;

%% Do optimisation - repeat frame-by-frame analysis
for n=1:length(GCI)
    
    % get framing information
    if n==1
        pulseLen=round((GCI(n+1)-GCI(n))*pulseNum);
        F0_cur=fs/(round(GCI(n+1)-GCI(n)));
    else pulseLen=round((GCI(n)-GCI(n-1))*pulseNum);
        F0_cur=fs/(round(GCI(n)-GCI(n-1)));
    end
    
     if GCI(n)-round(pulseLen/2) > 0
        start=GCI(n)-round(pulseLen/2);
    else start=1;
    end
    
    if GCI(n)+round(pulseLen/2) <= length(glot)
        finish = GCI(n)+round(pulseLen/2);
    else finish = length(glot);
    end
    
    % Get fft spectrum
    glot_seg=glot(start:finish).*hanning(finish-start+1);
    glot_seg=glot_seg(:);
    glot_seg_spec=20*log10(abs(fft(glot_seg)));
    
    % Get predicted R parameters from Rd
    [Ra_cur,Rk_cur,Rg_cur] = Rd2R(Rd_opt(n),EE(n),F0_cur);
    Ra_cur=0.01;

    if sum(glot_seg)==0
        params.EE(n)=EE(n);
        params.F0(n)=F0_cur;
        params.Ra(n)=0;
        params.Rk(n)=0;
        params.Rg(n)=0;
        
    else
        optFcn = @(R) get_time_err(glot_seg,R(1),R(2),R(3),EE(n),F0_cur, ...
            fs,GCI(n)-start,finish-GCI(n),pulseLen,glot_seg_spec,MVF);
        try
            R = fminsearch(optFcn,[Ra_cur Rk_cur Rg_cur],options);
        catch R(1)=Ra_cur; R(2)=Rk_cur; R(3)=Rg_cur;
        end

        params.EE(n)=EE(n);
        params.F0(n)=F0_cur;
        params.Ra(n)=R(1);
        params.Rk(n)=R(2);
        params.Rg(n)=R(3);
        
        pulse = lf_cont(params.F0(n),fs,params.Ra(n),params.Rk(n),params.Rg(n),params.EE(n));
        pulse_int=integrat(pulse,fs);
        params.UP(n)=max(pulse_int);
    end
end

params.OQ=(1+params.Rk)./(2*params.Rg);
params.Rd=params.Rd(:)';

%% Do plots 
t=(1:length(glot))/fs;
if nargin >3
    if plots
        figure
        subplot(211), 
        stem(t(GCI),glot(GCI),'k'), legend('GCIs'), hold on
        plot(t,glot), ylim([-.3 .3]), xlim([t(1) t(length(glot))]),
        hold on, stem(t(GCI),glot(GCI),'k')
        title('Estimated voice source signal','FontSize',16,'FontWeight','bold')
        xlabel('Time (seconds)','FontSize',14)
        ylabel('Amplitude','FontSize',14)
        set(gca,'FontSize',12)
        subplot(212), plot(t(GCI),Rd_n','rx'), xlim([t(1) t(length(glot))]),
        hold on, 
        for n=1:length(GCI)-1
            plot([t(GCI(n)) t(GCI(n+1))],[Rd_opt(n) Rd_opt(n+1)],'linewidth',2)
        end
        title('Optimal path of Rd values','FontSize',16,'FontWeight','bold')
        xlabel('Time (seconds)','FontSize',14)
        ylabel('Rd','FontSize',14)
        set(gca,'FontSize',12)
       % subplot(313), plot(GCI,params.OQ),xlim([0 length(glot)]), ylim([0 1])
    end
end

function err = get_time_err(glot_seg,Ra,Rk,Rg,EE,F0,fs,start,finish,pulseLen,glot_seg_spec,MVF)

if isinf(Rg)
    Rg=1;
end
pulse = lf_cont(F0,fs,abs(Ra),abs(Rk),abs(Rg),EE);
LFgroup = makePulseCentGCI(pulse,pulseLen,start,finish);
LFgroup_win=LFgroup(:).*hanning(length(LFgroup));
LFgroup_win=LFgroup_win(:);
freq=linspace(0,fs,length(glot_seg));
LFgroup_win_spec=20*log10(abs(fft(LFgroup_win)));

% Time domain error function
cor_time = corrcoef(glot_seg,LFgroup_win);
if isnan(cor_time)
    err_time = Inf;
else
    cor_time=abs(cor_time(2));
    err_time=1-cor_time;
end

        
% Frequency domain error function
cor_freq = corrcoef(glot_seg_spec(freq<MVF),LFgroup_win_spec(freq<MVF));
cor_freq=abs(cor_freq(2));
err_freq=1-cor_freq;

err=err_time+err_freq;