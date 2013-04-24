function GCI = GCI_creak_postproc(GCI,creak,search_reg,rep,removeThresh,repNum)

% Function to do the post-processing step to remove false positive GCIs
% detected in creaky voice regions

% Separate GCIs detected in creaky voice regions from other regions
creak_GCI=GCI(creak(GCI)==1);
GCI(creak(GCI)==1)=[];
    
for m=1:repNum
    n=2;
    while n < length(creak_GCI)
        cur_rep_max = abs(min(rep(creak_GCI(n)-round(search_reg):creak_GCI(n)+round(search_reg))));

        if mean([abs(rep(creak_GCI(n-1))) abs(rep(creak_GCI(n+1)))])*removeThresh > cur_rep_max
             creak_GCI(n)=NaN;
            n=n+2;
        else n=n+1;
        end
    end
    creak_GCI(isnan(creak_GCI))=[];
end
    
GCI=sort(unique([GCI(:)' creak_GCI(:)']));