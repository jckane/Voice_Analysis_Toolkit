function [sig,start,finish,UP] = genLFsig(params,glot,GCI,fs)

% Function to generate synthetic source signal using LF model based
% parameter values.
%
% USAGE:    
%       Input:
%             params  : struct containing LF model based parameters
%             glot   : glottal flow derivative estimate
%             GCI : Glottal closure instants (in samples)
%             fs : sampling frequency (in Hz)
%
%
%       Output:
%             y  : filtered signal
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = LFparamCheck(params);

sig = zeros(1,length(glot));
start = zeros(1,length(GCI));
finish = zeros(1,length(GCI));
UP = zeros(1,length(GCI));
maxCnt=10;

F0min=20;
F0max=500;

%% Generate synthetic source signal
for n=1:length(GCI)
        if params.F0(n) > F0min && params.F0(n) < F0max
            %[params.F0(n),fs,params.Ra(n),params.Rk(n),params.Rg(n),params.EE(n)]
            pulse = lf_cont(params.F0(n),fs,params.Ra(n),params.Rk(n),params.Rg(n),params.EE(n));
            pulse_int=integrat(pulse,fs);
            UP(n)=max(pulse_int);
            cnt=1;
            while any(isnan(pulse)) && cnt < maxCnt
                params.Rg(n)=params.Rg(n)+0.01;
                pulse = lf_cont(params.F0(n),fs,params.Ra(n),params.Rk(n),params.Rg(n),params.EE(n));
                cnt=cnt+1;
            end
            if cnt==maxCnt
                pulse=zeros(1,length(pulse));
            end
            [~,idx]=min(pulse);
            start(n)=GCI(n)-idx-1;
            finish(n)=start(n)+length(pulse)-1;
            if start(n) > 0 && finish(n) <=length(sig)
                sig(start(n):finish(n)) = sig(start(n):finish(n))+pulse;
            end
        end
end

