function [T1,T2] = findAmid_t(glot_adj,Amid,Tz)

% Function to find the start and stop positions of the quasi-open phase.

T1=0;
T2=0;
if Tz~=0
    n=Tz;

    while glot_adj(n) > Amid && n > 3
        n=n-1;
    end
    T1=n;
    n=Tz;

    while glot_adj(n) > Amid && n < length(glot_adj)-2
        n=n+1;
    end
    T2=n;
end