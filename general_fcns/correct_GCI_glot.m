function GCI_out = correct_GCI_glot(glot,GCI,fs)

% Function to adjust GCI locations to coincide with glottal source
% excitation minima
%
% USAGE:    
%       Input:
%             glot   : glottal flow derivative estimate
%             GCI : Glottal closure instants (in samples)
%
%       Output:
%             GCI_out: adjusted GCIs
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initial settings
search_rate=0.2;
GCI_out=zeros(1,length(GCI));
maxPulseLen=.03;

%% Do processing
if length(GCI)>1
    for n=1:length(GCI)

        if n==1
            pulseLen=GCI(n+1)-GCI(n);
        else pulseLen=GCI(n)-GCI(n-1);
        end
        
        if pulseLen < maxPulseLen*fs
            if GCI(n) - round(pulseLen*search_rate)>0
                start = GCI(n) - round(pulseLen*search_rate);
            else start=1;
            end
            if GCI(n) + round(pulseLen*search_rate)<=length(glot)
                finish = GCI(n) + round(pulseLen*search_rate);
            else finish = length(glot);
            end

            [~,idx]=min(glot(start:finish));
            if isempty(idx)==0
                GCI_out(n)=idx+start-1;
            else GCI_out(n)=NaN;
            end
        else GCI_out(n)=GCI(n);
        end
        

    end
end
GCI_out(isnan(GCI_out))=[];

    