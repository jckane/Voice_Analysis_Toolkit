function [f0_inter,f0_samp] = create_continuous_smooth_f0(f0,VUV,x)

% Function to create and f0 contour which is interpolated over unvoiced
% regions and is heavily smoothed with median and moving average filters.

%% Initial settings
F0_min=50;
F0_max=500;
med_len=17;
sm_len=17;
F0min=min(medfilt1(f0(VUV==1&f0>F0_min&f0<F0_max),med_len));
F0max=max(medfilt1(f0(VUV==1&f0>F0_min&f0<F0_max),med_len));

f0(f0<F0min)=F0min;
f0(f0>F0max)=F0max;

VUV2=zeros(1,length(VUV));
VUV2(VUV==0)=1;
VUV=VUV2;

VUV_reg=bin2reg(VUV,VUV);

VUV_reg(1,1)=1;
N=size(VUV_reg,1);
f0_inter=f0;

f0_inter(1:VUV_reg(1,2)+1)=f0_inter(VUV_reg(1,2)+10);

for n=2:N
   
    start=VUV_reg(n,1)-1;
    stop=VUV_reg(n,2)+1;
    f0_int_cur=interp1([start stop],f0([start stop]),start:stop);
    f0_inter(start:stop)=f0_int_cur;
end

f0_inter=smooth(medfilt1(f0_inter,med_len),sm_len);
f0_samp=interp1(linspace(1,length(x),length(f0_inter)),f0_inter,1:length(x));
