function [Ra,Rk,Rg] = Rd2R(Rd,EE,F0)

% Function derive default R parameter settings from a given Rd value.
%
% USAGE:    
%       Input:
%             Rd : Rd value
%             EE : EE (excitation strength)
%             f0 : fundamental frequency (in Hz)
%
%       Output:
%             Ra :
%             Rk :
%             Rg :
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ra = (-1+(4.8*Rd))/100;
Rk = (22.4 +(11.8*Rd))/100;
EI = (pi*Rk*EE)/2;
UP = (Rd*EE)/(10*F0);
Rg = EI/(F0*UP*pi);