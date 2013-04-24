function params = create_LFparams_zeros(gci)

% Function to create dummy LF model parameter arrays

N=length(gci);
params.Ra=zeros(1,N);
params.Rk=zeros(1,N);
params.Rg=zeros(1,N);
params.Rd=zeros(1,N);
params.EE=zeros(1,N);
params.UP=zeros(1,N);
params.OQ=zeros(1,N);
params.Fa=zeros(1,N);
params.F0=zeros(1,N);