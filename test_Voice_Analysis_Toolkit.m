function test_Voice_Analysis_Toolkit(x,fs)

% Function to apply the algorithms of the Voice Analysis Toolkit to an
% inputted speech signal and plot the results.

addpath(genpath(pwd));

%% Initial settings
F0min=20;
F0max=500;
med_len=5;
sm_len=5;
creak_pp_thresh=0.3;

% Plotting settings
fig_left=250;
fig_bottom=500;
fig_width=1300;
fig_height=300;
y_min=-1;
y_max=1;
y2_min=0;
y2_max=.15;
tick_step_num=5;

%% f0 and VUV detection using SRH algorithm (from Drugman's GLOAT toolkit)
disp('f0 and voicing detection using SRH algorithm (GLOAT toolkit)')
[f0,VUV] = SRH_PitchTracking(x,fs,F0min,F0max);
F0med=median(f0(f0>F0min&f0<F0max&VUV==1));
VUV_inter=interp1(linspace(1,length(x),length(VUV)),VUV,1:length(x));
VUV_inter(VUV_inter<0.5)=0;
VUV_inter(VUV_inter>=0.5)=1;

%% Creaky voice detection
disp('Detecting creaky voice using Kane-Drugman method')
[creak_pp,creak_bin,creak_t] = CreakyDetection_CompleteDetection(x,fs);
creak_pp_inter=interp1(linspace(1,length(x),length(creak_pp)),creak_pp,1:length(x));
creak_bin_inter=interp1(linspace(1,length(x),length(creak_bin)),creak_bin,1:length(x));
creak_bin_inter(creak_bin_inter<.5)=0;creak_bin_inter(creak_bin_inter>=.5)=1;
VUV_inter(creak_pp_inter>creak_pp_thresh)=1;

%% GCI detection - SQ-VQ algorithm
disp('GCI detection using SE-VQ algorithm')
%[GCI,~,res] = SE_VQ(x,fs,f0,VUV,creak_bin_inter); % Standard SE-VQ algorithm
[GCI,~,res] = SE_VQ_varF0(x,fs,f0,VUV,creak_bin_inter); % SE-VQ for highly varying f0
GCI=unique(sort(GCI));
GCI(VUV_inter(GCI)==0)=[];

%% Glottal inverse filtering using IAIF algorithm (Alku et al. 1992)
disp('Glottal inverse filtering using IAIF algorithm (Alku et al. 1992)')
g_iaif=IAIF(x,fs,GCI);
g_iaif=g_iaif/max(abs(g_iaif));

%% Derive peakSlope and glottal-peakSlope
peakSlope = get_peakSlope(x,fs);
peakSlope_g = get_peakSlope(g_iaif,fs);

%% Extracting MDQ parameter
disp('Extracting MDQ parameter')
MDQ = get_MDQ(res,fs,GCI);
MDQ_frame = get_MDQ_frame(res,fs,F0med);
MDQ=smooth(medfilt1(MDQ,med_len),sm_len);
MDQ_frame=smooth(medfilt1(MDQ_frame,med_len),sm_len);

MDQ_inter=interp1(linspace(1,length(x),length(MDQ)),MDQ,1:length(x));
MDQ_inter(VUV_inter==0)=0;

%% Do plots
close
figure('Position',[fig_left, fig_bottom, fig_width, fig_height])
t=(1:length(x))/fs;
g_iaif=g_iaif/max(abs(g_iaif));
[AX,H1,H2]=plotyy(t,g_iaif,t,MDQ_inter);

set(AX(1),'YLim',[y_min y_max])
set(AX(2),'YLim',[y2_min y2_max])
set(AX(1),'YTick',linspace(y_min,y_max,tick_step_num))
set(AX(2),'YTick',linspace(y2_min,y2_max,tick_step_num))
set(get(AX(1),'Ylabel'),'String','Amplitude')
set(get(AX(2),'Ylabel'),'String','MDQ')
xlabel('Time (seconds)','FontSize',14)
hold on,
stem(t(GCI),ones(1,length(GCI))*-.3,'k')
plot(creak_t,creak_pp,'--m','LineWidth',2)
legend('Speech waveform','GCIs','Creaky voice probability')
drawnow
    
