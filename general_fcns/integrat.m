function y = integrat(x,Fs)

% [y,H,P] = integrat(Fs,x)
% Function to apply a simple integrator

y=zeros(1,length(x));

Ts=1/Fs;
y(1)=Ts*x(1);

% Difference equation
for n=2:length(x),
    y(n)=(Ts*x(n))+y(n-1);
end