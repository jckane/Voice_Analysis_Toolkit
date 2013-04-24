function LFgroup = makePulseCentGCI(pulse,winLen,start,finish)

% Function to get a GCI centred LF pulse group
%
% USAGE:    
%       Input:
%             pulse   : LF model pulse
%             winLen : length of the window
%
%       Output:
%             LFgroup: GCI centred group of LF pulses
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pulse=pulse(:)';
[~,idx]=min(pulse);
pulseLen = length(pulse);
group_idx=idx+pulseLen;
pulseGroup=[pulse pulse pulse];

if nargin < 3
    if rem(winLen,2)~=0
        start = group_idx-ceil(winLen/2);
    else start = group_idx-winLen/2;
    end
    finish = group_idx+floor(winLen/2);
else start = group_idx-start;
    finish = group_idx+finish;
end

if finish > length(pulseGroup) || start < 1
    LFgroup = zeros(1,finish-start+1);
else LFgroup=pulseGroup(start:finish);
end