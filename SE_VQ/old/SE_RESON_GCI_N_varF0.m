function [gci,goi,Average,gci_relAmp,Interval] = ...
    SE_RESON_GCI_N_varF0(signal,Fs,res,T0mean,Ncand)

% Thomas Drugmans SEDREAMS method for detecting glottal closure and glottal
% opening instants using information from the mean based signal to
% determine a search area and then peak picking the from the LPC residual
%
% Method is slightly modified in order to select the N strongest
% LP-residual peaks in each interval. Peaks are normalised by the maximum
% peak and subtracted from 1 to be used a target cost in a dynamic
% programming method.

% USAGE: INPUT
%        signal - speech signal
%        Fs     - sampling frequency
%        res    - LPC residual, obtained using GetLPCresidual function
%
%        OUTPUT
%        gci    - glottal closure instant locations in samples
%        goi    - glottal opening instant locations in samples


%%%%%%%%%%%%%%%
gci=[];
goi=[];
gci_relAmp=[];

if nargin < 5
    Ncand=5; % Number of candidate GCIs to detect in each interval
end

% Obtain the mean-based signal
Average=zeros(1,length(signal));
Interval=zeros(1,length(signal));

halfL=round((1.7*T0mean(1))/2);
BlackWin=blackman(2*halfL+1);

Maxis=[];
Minis=[];
PosAverage=zeros(1,length(signal));
StepExp=4;
Step=2^StepExp;

%% Calculate mean-based signal
for m=halfL+1:Step:length(signal)-halfL
    
    halfL=round((1.7*T0mean(m))/2);
    BlackWin=blackman(2*halfL+1);
    
    start=m-halfL;
    stop=m+halfL;
    if stop > length(signal)
        break
    end
    
    vec=signal(m-halfL:m+halfL);
    
    vec=vec.*BlackWin;
    Average(m)=mean(vec);
    PosAverage(m)=1;
end

posis=find(PosAverage==1);
for tmp=2:length(posis)-1
    if (Average(posis(tmp))>Average(posis(tmp-1))) && (Average(posis(tmp))>Average(posis(tmp+1)))
        Maxis=[Maxis posis(tmp)];
    end
    
    if (Average(posis(tmp))<Average(posis(tmp-1))) && (Average(posis(tmp))<Average(posis(tmp+1)))
        Minis=[Minis posis(tmp)];
    end
end


for StepExp=3:-1:0
    Step=2^StepExp;
    
    for tmp=1:length(Maxis)
        m=Maxis(tmp)-Step;
        
         halfL=round((1.7*T0mean(m))/2);
         BlackWin=blackman(2*halfL+1);
         
         start=m-halfL;
        stop=m+halfL;
        if stop > length(signal)
            break
        end
        
        vec=signal(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        Average(m)=mean(vec);
        PosAverage(m)=1;
        
        m=Maxis(tmp)+Step;
        vec=signal(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        Average(m)=mean(vec);
        PosAverage(m)=1;
        
        VecTmp=[Average(Maxis(tmp)-Step) Average(Maxis(tmp)) Average(Maxis(tmp)+Step)];
        [~,PosTmp]=max(VecTmp);
        Maxis(tmp)=Maxis(tmp)+(PosTmp-2)*Step;
    end
    
    for tmp=1:length(Minis)
        m=Minis(tmp)-Step;
        
        halfL=round((1.7*T0mean(m))/2);
         BlackWin=blackman(2*halfL+1);
         
         start=m-halfL;
        stop=m+halfL;
        if stop > length(signal)
            break
        end
        
        vec=signal(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        Average(m)=mean(vec);
        PosAverage(m)=1;
        
        m=Minis(tmp)+Step;
        vec=signal(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        Average(m)=mean(vec);
        PosAverage(m)=1;
        
        VecTmp=[Average(Minis(tmp)-Step) Average(Minis(tmp)) Average(Minis(tmp)+Step)];
        [~,PosTmp]=min(VecTmp);
        Minis(tmp)=Minis(tmp)+(PosTmp-2)*Step;
    end
    
end


Average=Average/max(abs(Average));

%%
if isempty(Maxis) || isempty(Minis) || length(Maxis) <2 || length(Minis) <2 
    return
end
if Maxis(1)<Minis(1)
    Maxis(1)=[];
end


if Minis(end)>Maxis(end)
    Minis(end)=[];
end


%% Determine the median position of GCIs within the cycle
res=res/max(abs(res));
%Posis=find((res)>0.4);
Posis=find(abs(res)>0.4);
RelPosis=zeros(length(Posis),1);
for k=1:length(Posis)

    Dists=abs(Minis-Posis(k));
    [~,pos]=min(Dists);

    interv=Maxis(pos)-Minis(pos);
    RelPosis(k)=(Posis(k)-Minis(pos))/interv;
end


RatioGCI=median(RelPosis);


% Detect GCIs and GOIs
gci=zeros(Ncand,length(Minis));
gci_relAmp=zeros(Ncand,length(Minis));
goi=zeros(1,length(Minis));

for k=1:length(Minis)
    interv=Maxis(k)-Minis(k);
    %     if type(Minis(k)+round(interv/2))~=0
    alpha=RatioGCI-0.25;
    start=Minis(k)+round(alpha*interv);
    alpha=RatioGCI+0.35;
    stop=Minis(k)+round(alpha*interv);

    if start<1
        start=1;
    end
    if stop>length(res)
        stop=length(res);
    end
    Interval(start:stop)=1;
    vec=res(start:stop);
    if stop < start
        break
    end
    % Get n GCI candidates within defined interval
    maxi=zeros(Ncand,1);
    posi=zeros(Ncand,1);
    for m=1:Ncand
        [maxi(m),posi(m)]=max(abs(vec));
        vec(posi(m))=0;
        gci(m,k)=start+posi(m)-1;
        %gci_relAmp(m,k)=maxi(1)/maxi(m);
        gci_relAmp(m,k)=1 - (maxi(m)/maxi(1));
    end

    alpha=RatioGCI+0.4;
    start=Minis(k)+round(alpha*interv);
    alpha=RatioGCI+1.5;
    stop=Minis(k)+round(alpha*interv);

    if start<1
        start=1;
    end
    if stop>length(res)
        stop=length(res);
    end

    vec=res(start:stop);
    [~,posi]=max(abs(vec));
    goi=[goi start+posi-1];
    %     end
end
