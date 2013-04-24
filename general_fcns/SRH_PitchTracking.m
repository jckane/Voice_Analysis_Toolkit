function [f0,VUVDecisions,SRHVal,time] = SRH_PitchTracking(wave,Fs,F0min,F0max)

% USAGE:
% [f0,VUVDecisions,SRHVal] = SRH_PitchTracking(wave,Fs,F0min,F0max)
% 
% INPUTS:
%     - wave is the speech signal
%     - Fs is the sampling frequency (Hz)
%     - F0min is the minimum value for F0 search (Hz)
%     - F0max is the maximum value for F0 search (Hz)
%    
% OUPUTS:
%     - f0 is the vector of F0 values (with an hopsize of 10ms)
%     - VUVDecisions is the vector of voiced-unvoiced VUVDecisions (with an hopsize of 10ms)
%     - SRHVal is the vector of SRH values (with an hopsize of 10ms)
%     
% For more information, please read:
% 
% T.Drugman, A.Alwan, "Joint Robust Voicing Detection and Pitch Estimation
% Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011
%
% Please refer to this work in your publication if you use this code.
% 
% Code written by Thomas Drugman in TCTS Lab, University of Mons, Belgium.
% 
% Copyright (C) 2000-2011 Thomas Drugman - TCTS Lab
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Matlab version used for developing this code: 7.0
%
% http://tcts.fpms.ac.be/~drugman/      thomas.drugman@umons.ac.be

%display('Tracking Pitch with SRH Method')
pause(0.0001)

if Fs>16000
    wave=resample(wave,16000,Fs);
    Fs=16000;
end

LPCorder=round(3/4*Fs/1000);
Niter=2;

res = GetLPCresidual(wave,round(25/1000*Fs),round(5/1000*Fs),LPCorder);

%% Estimate the pitch track in 2 iterations
% PrevF0meanEst=0;
for Iter=1:Niter   

    [f0,SRHVal,time] = SRH_EstimatePitch(res',Fs,F0min,F0max);
    
    posiTmp=find(SRHVal>0.1);
    
    if length(posiTmp)>1
        F0meanEst=median(f0(posiTmp));
        
        % Delta=abs(F0meanEst-PrevF0meanEst);
        F0min=round(0.5*F0meanEst);
        F0max=round(2*F0meanEst);
        
        % PrevF0meanEst=F0meanEst;
    end
    
end

% Voiced-Unvoiced decisions are derived from the value of SRH (Summation of
% Residual Harmonics)
VUVDecisions=zeros(1,length(f0));
Threshold=0.07;
if std(SRHVal)>0.05
    Threshold=0.085;
end

pos= SRHVal>Threshold;
VUVDecisions(pos)=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f0,SRHVal,time] = SRH_EstimatePitch(sig,Fs,F0min,F0max)

time=[];
start=1;
stop=round(100/1000*Fs);
shift=round(10/1000*Fs);

Nframes=floor((length(sig)-stop)/shift) + 1;
SRHTot=zeros(Nframes,F0max);
F0s = zeros(1,Nframes);
SRHVal = zeros(1,Nframes);

BlackWin=blackman(stop-start+1);

index=1;
while stop<=length(sig)
    time(index) = start;
    seg=sig(start:stop);
    seg=seg.*BlackWin;    
    seg=seg-mean(seg);

    Spec=fft(seg,Fs);
    Spec=abs(Spec(1:Fs/2));    
    Spec=Spec/sqrt(sum(Spec.^2));
        
    SRHs=zeros(1,F0max);    
    
    % SRH spectral criterion
    for freq=F0min:F0max
        SRHs(freq)=(Spec(freq)+Spec(2*freq)+Spec(3*freq)+Spec(4*freq)+Spec(5*freq))-(Spec(round(1.5*freq))+Spec(round(2.5*freq))+Spec(round(3.5*freq))+Spec(round(4.5*freq)));
    end
    
    SRHTot(index,:)=SRHs;

    [maxi,posi]=max(SRHs);
    F0frame=posi;
        
    F0s(index)=F0frame;
    SRHVal(index)=SRHs(F0frame);
    
    start=start+shift;
    stop=stop+shift;
    index=index+1;
end

f0=[0 0 0 0 0 F0s 0 0 0 0 0];
SRHVal=[0 0 0 0 0 SRHVal 0 0 0 0 0];
