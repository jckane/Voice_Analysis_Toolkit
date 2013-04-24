function y = RCVD_reson_GCI(res,fs,F0mean)

% Function to use the resonator used in RCVD (creaky voice detection), 
% applied to the LP-residual signal and give output


%% Configure resonator (using settings in RCVD)
Phi=2*pi*1*F0mean/fs;
Rho=0.9; % Set to a narrow bandwidth
rep=filtfilt([1 0 0],[1 -2*Rho*cos(Phi) Rho^2],res); % Filter forwards and backwards to have zero-phase
y=rep/max(abs(rep));