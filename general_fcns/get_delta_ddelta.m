function [X_d,X_dd] = get_delta_ddelta(X)

% Function to get delta and delta-delta coefficients of input feature
% matrix.

X_d=zeros(size(X));
X_dd=zeros(size(X));

% delta
for n=1:size(X,2)
    X_d(2:end,n)=diff(X(:,n));
end

% delta-delta
for n=1:size(X,2)
    X_dd(2:end,n)=diff(X_d(:,n));
end
