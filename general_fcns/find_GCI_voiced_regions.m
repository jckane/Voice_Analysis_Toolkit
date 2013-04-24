function [VUV_reg,GCI] = find_GCI_voiced_regions(VUV,GCI,minGCI)

% Function to find voicing regions containing more than minGCI GCI samples

%% Initial settings
VUV_regions = bin2reg(VUV,VUV);
GCI(VUV(GCI)==0)=[];
N=size(VUV_regions,1);
VUV_reg=[];

if nargin < 3
    minGCI=3;
end

%% Do processing
for n=1:N
   
    numGCI=length(GCI(GCI>=VUV_regions(n,1)&GCI<=VUV_regions(n,2)));
    
    if numGCI>=minGCI
        VUV_reg=[VUV_reg; VUV_regions(n,:)];
    end
end

