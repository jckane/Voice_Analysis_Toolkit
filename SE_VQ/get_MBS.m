function MBS = get_MBS(x,fs,T0mean)

% Function to calculate the mean-based signal as is done in SEDREAMs

% Obtain the mean-based signal
MBS=zeros(1,length(x));

halfL=round((1.6*T0mean(1))/2);

StepExp=3;
Step=2^StepExp;

%% Calculate mean-based signal
for m=halfL+1:Step:length(x)-halfL
    
    if length(T0mean)==1
        halfL=round((1.7*T0mean(1))/2);
    else halfL=round((1.7*T0mean(m))/2);
    end
    BlackWin=blackman(2*halfL+1);
    
    start=round(m-halfL);
    stop=round(m+halfL);
    if stop > length(x)
        break
    end

    if start > 0
        vec=x(start:stop);

        vec=vec.*BlackWin;
        MBS(m)=mean(vec);
    end
end

t=find(MBS~=0);
MBS=interp1(t,MBS(t),1:length(x));
MBS(isnan(MBS))=0;
MBS=zeroPhaseHPFilt(MBS,fs,70,10,0);
MBS=MBS/max(MBS);
MBS=smooth(MBS,7);

