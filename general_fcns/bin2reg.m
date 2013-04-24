function regions = bin2reg(regions_bin,x)

% Function to convert a binary region vector to a vector of separate start
% and end points

% USAGE:    
%       Input:
%             regions_bin : binary region vector
%             x  : input signal
%
%       Output:
%             regions : 2D vector with start and end times for each region.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

regions=[];

n=1;
m=0;

while n < length(regions_bin)
    if regions_bin(n)==1 
        m=m+1;
        regions(m,1) = n;
        
        while regions_bin(n)==1 && n <length(x)
            n=n+1;
        end
        
        regions(m,2) = n-1;
    else n=n+1;
    end

end