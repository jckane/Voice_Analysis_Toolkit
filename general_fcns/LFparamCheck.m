function params = LFparamCheck(params)

% Function to check parameters for spurious values and correcting them

lb = [0.001 0.05 0.65]; % parameter value lower bounds
ub = [0.2 0.5 2]; % parameter value upper bounds
    

% Do processing
for n=1:length(params.EE)
    
    % Do checks to ensure sensible parameter values
    if params.Ra(n) < lb(1)
        params.Ra(n) = lb(1);
    elseif params.Ra(n) > ub(1)
        params.Ra(n) = ub(1);
    end
    
    if params.Rk(n) < lb(2)
        params.Rk(n) = lb(2);
    elseif params.Rk(n) > ub(2)
        params.Rk(n) = ub(2);
    end
    
    if params.Rg(n) < lb(3)
        params.Rg(n) = lb(3);
    elseif params.Rg(n) > ub(3)
        params.Rg(n) = ub(3);
    end
    
end