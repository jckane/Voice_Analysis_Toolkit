function [H2H1,f0] = get_creak_H2H1(res,fs,F0mean)

% Function to calculate H2-H1 contour based on resonator applied to
% LP-residual
%
% REFERENCE:
%       Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%
%       Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', 27(4), pp. 1028-1047, Computer Speech and Language.
%
%       Drugman, T., Kane, J., Gobl, C., (2012) `Resonator-based creaky 
%       voice detection', Interspeech 2012, Portland, Oregon, USA.

%% Initial settings
buffer=0;
Start=1;
Stop=50/1000*fs;

Hop=10/1000*fs;
Win=hanning(Stop);
Ind=1;
f0=[];
Delta=[];

moving_average_size=100/1000*fs;
moving_average_frame=moving_average_size/Hop;

% Resonator 2 settings
Phi=2*pi*1*F0mean/fs;
Rho=0.97;
rep=filter([1 0 0],[1 -2*Rho*cos(Phi) Rho^2],res);

rep=rep/max(abs(rep));

% Resonator 1 settings
Phi=2*pi*1*F0mean/fs;
Rho=0.8;
rep2=filter([1 0 0],[1 -2*Rho*cos(Phi) Rho^2],res);

%% Do processing
while Stop<length(rep)
    Sig=rep2(Start:Stop)';    
    Sig=Sig.*Win;
    
    Corrs=xcorr(Sig,Sig);
    Corrs(1:length(Sig))=[];
    
    for k=1:length(Corrs)
       Corrs(k)=Corrs(k)*(length(Corrs)/(length(Corrs)-(k-1))); 
    end

    Corrs2=Corrs;
    k=1;
    while (Corrs2(k)>0)
        Corrs2(k)=0;
        k=k+1;
        if k==length(Corrs2)
            break;
        end
    end
    [maxi,posi]=max(Corrs2);
    
    
    if posi<5
        posi=5;
    end    
    
    F0=round(fs/posi);
    
    f0(Ind)=F0;
    F0=round(f0(Ind));
    
    Sig=rep(Start:Stop)';
    Sig=Sig.*Win;
    
    Corrs=xcorr(Sig,Sig);
    Corrs(1:length(Sig))=[]; Sig=Sig.*Win;
    
    Corrs=xcorr(Sig,Sig);
    Corrs(1:length(Sig))=[];
    
    
    
    Spec=fft(Corrs,fs);
%     Spec=fft(Sig,fs);
    Spec=abs(Spec(1:fs/2));
    
    Spec=Spec/sqrt(sum(Spec.^2));
    Spec=20*log10(Spec);
        
    Delta(Ind)=(Spec(2*F0)-Spec(F0))+buffer;

    
    Start=Start+Hop;
    Stop=Stop+Hop;
    Ind=Ind+1;
end

%% Interpolate to update every sample
if exist('Delta') 
    if isempty(Delta)==0
        Delta=smooth(Delta,moving_average_frame);
        f0=medfilt1(f0,7);
        f0=smooth(f0,moving_average_frame);
        Delta=interp1(linspace(1,length(res),length(Delta)),Delta,1:length(res));
        %Delta=resample(Delta,length(res),length(Delta));
        Delta=Delta-buffer;
        
        f0=interp1(linspace(1,length(res),length(f0)),f0,1:length(res));

        %Delta2=smooth(Delta,moving_average_size); % 100ms moving average filter
        %f0=smooth(f0,moving_average_size); % 100ms moving average filter
        Delta2=Delta;

        H2H1=Delta2;
    else H2H1=zeros(1,length(res))-buffer;
        f0=zeros(1,length(res));
    end
else H2H1=zeros(1,length(res))-buffer;
    f0=zeros(1,length(res));
end