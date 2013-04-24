function y = nearest_neighbour_resample(x,x_idx,fs,idx,interval)
    
% Function to resample a parameter trajectory using a nearest neighbour
% approach. This is particularly useful for GCI synchronous parameters
% Input parameters:
%       x_idx = usually is the GCI locations
%       interval =  the maximum duration for searching around the location,
%                   otherwise output 0.

if nargin < 5
    interval=30/1000*fs;
end

if length(x_idx) > length(x)
    x_idx=x_idx(1:length(x));
    disp('Difference in nearest neighbour search!!')
end

y=zeros(1,length(idx));

for n=1:length(idx)
   
    [~,closest_idx]=min(abs(x_idx-idx(n)));
    dist_closest = abs(x_idx(closest_idx)-idx(n));
    
    if dist_closest < interval
        y(n) = x(closest_idx);
    end
end
