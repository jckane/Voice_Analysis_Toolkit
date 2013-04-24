function ss = get_spec_stat(x,fs,GCI)

% Function to get spectral stationarity measure using the same method as
% used by Talkin in fxrapt/

med_len=7;
sm_len=7;

sp=filter([1 -.99],1,x); % preemphasised speech for LPC calculation
lpcord=2+round(fs/1000);

frameLen = 25/1000*fs;
frameShift = 10/1000*fs;

ss=zeros(round((length(x)-frameLen)/frameShift),1);

start=1;
stop=start+frameLen-1;
seg_prev = sp(start:stop);
seg_prev_win = seg_prev(:).*hamming(length(seg_prev));

%% Do procesing
cnt=1;
while stop <= length(x)
    
    % get framing info
    start=start+frameShift-1;
    stop=start+frameLen-1;
    
    if stop > length(sp)
        break
    end
    
    seg_cur = sp(start:stop);
    seg_cur_win = seg_cur(:).*hamming(length(seg_cur));
    
    % Calculate spectral stationarity
    ss(cnt)=0.2/(distitar(lpcauto(seg_cur_win,lpcord),lpcauto(seg_prev_win,lpcord),'e')-0.8);

    seg_prev_win=seg_cur_win;
    cnt=cnt+1;
end

ss=smooth(medfilt1(ss,med_len),sm_len);
ss_inter=interp1(linspace(1,length(x),length(ss)),ss,1:length(x));
ss=ss_inter(GCI);