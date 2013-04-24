function [Outs,Decs] = CreakyDetection_DoClassification(FeatTot,ANN)

% Function to use the extracted features as part of the creak
% ANN classification
%
% REFERENCE:
%       Drugman, T., Kane, J., Gobl, C., `Automatic Analysis of Creaky
%       Excitation Patterns', Submitted to Computer Speech and
%       Language.
%
%       Kane, J., Drugman, T., Gobl, C., (2013) `Improved automatic 
%       detection of creak', 27(4), pp. 1028-1047, Computer Speech and Language.
%
%       Drugman, T., Kane, J., Gobl, C., (2012) `Resonator-based creaky 
%       voice det

%% Initial settings
Maxis=ANN.Maxis.Maxis;
Minis=ANN.Minis.Minis;
net=ANN.net.net;

ANN_decision_threshold=0.3;
med_len=3; % Length of median filtering

%% Do Z-score normalisation
m=size(FeatTot);
X=zeros(m(2),m(1));
for k=1:m(2)
    X(k,:)=FeatTot(:,k)';    
end

m=size(X);

for k=1:m(1)
    vec=X(k,:);
    
    mini=Minis(k);
    maxi=Maxis(k);
    
    pos=find(isnan(vec));
    vec(pos)=mini;
    
    X(k,:)=-1+((vec-mini)/(maxi-mini))*2;
end

%% Simulate ANN network
Outs=sim(net,X);

Outs2=medfilt1(Outs,med_len); % Median filtering to remove transients

Decs=zeros(1,length(Outs));
Decs(Outs2>ANN_decision_threshold)=1;
