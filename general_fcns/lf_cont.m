function g_LF = lf_cont(F0,fs,Ra,Rk,Rg,EE)

% Function to generate an LF model pulse given the Rparameter set.
% Optimisation algorithm for the solving the area balance is done using
% self-coded newton raphson method, so no in-built matlab functions used. 
% The C file lf_Area_newton.c needs to be mex-ed before this can be used.
% USAGE:    
%       Input:
%             F0   : glottal flow derivative estimate
%             fs  : sampling frequency (in Hz)
%             Ra : 
%             Rk : 
%             Rg : 
%             EE :
%
%       Output:
%             g_LF: LF model
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F0min=20;
F0max=500;

%% Set LF model parameters
T0=1/F0;
Ta=Ra*T0;
Te=((1+Rk)/(2*Rg))*T0; 
Tp=Te/(Rk+1);
Tb=((1-(Rk+1)/(2*Rg))*1/(F0));
Tc=Tb+Te; 

if F0<F0min || F0>F0max
    g_LF=NaN;
else


    %% Solve area balance using newton raphson method
    [alpha,epsi]=lf_Area_newton(Tc,fs,Tp,Te,Ta,EE);
    omega=pi/Tp;
    E0=-(abs(EE))/(exp(alpha*Te)*sin(omega*Te));

    %% Generate open phase and closed phase and combine
    T1=1/fs:1/fs:Te;
    T2=(T1(end)+1/fs:1/fs:Tc);
    g_op=E0*exp(alpha(end)*T1).*sin(omega*T1);
    g_cl=(-EE/(epsi*Ta))*(exp(-epsi*(T2-Te))-exp(-epsi*Tb));
    g_LF=[g_op g_cl];
end
