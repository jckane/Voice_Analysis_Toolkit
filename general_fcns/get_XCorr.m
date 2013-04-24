function Xcorr_amp=get_XCorr(x1,x2,fs)

% Function to cross-correlation between to signals on a fixed frame-basis.

frameLen=32/1000*fs;
frameShift=10/1000*fs;
buffer=2/1000*fs; % Exclusion regions for ACF peak search

Xcorr_amp=zeros(1,round((length(x1)-frameLen)/frameShift));

start=1;
stop=start+frameLen-1;

%% Do processing
cnt=1;
while stop < length(x1)
    
    % Get windowed frame
    frame1=x1(start:stop);
    frame_win1=frame1(:).*hamming(length(frame1));
    
    frame2=x2(start:stop);
    frame_win2=frame2(:).*hamming(length(frame2));
    
    % Get normalised autocorrelation function and peak above buffer
    [Xcorr,lags]=xcorr(frame_win1,frame_win2,'coeff');
%     subplot(211), plot(frame_win1), hold on, plot(frame_win2,'r'), hold off
%     subplot(212), plot(lags,Xcorr), 
    Xcorr_amp(cnt)=max(Xcorr);
%     Xcorr_amp(cnt)
%     pause(.5)
    
    % Get next frame info
    start=start+frameShift-1;
    stop=start+frameLen-1;
    cnt=cnt+1;
    
    
end
