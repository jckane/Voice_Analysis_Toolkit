function ACF_amp=get_ACF_amp(x,fs)

% Function to detect strongest peak in the autocorrelation function
% measured on a fixed frame-basis

frameLen=32/1000*fs;
frameShift=10/1000*fs;
buffer=2/1000*fs; % Exclusion regions for ACF peak search

ACF_amp=zeros(1,round((length(x)-frameLen)/frameShift));

start=1;
stop=start+frameLen-1;

%% Do processing
cnt=1;
while stop < length(x)
    
    % Get windowed frame
    frame=x(start:stop);
    frame_win=frame(:).*hamming(length(frame));
    
    % Get normalised autocorrelation function and peak above buffer
    [ACF,lags]=xcorr(frame_win,'coeff');
    ACF(lags<0)=0;
    ACF(lags<buffer)=0;
    ACF_amp(cnt)=max(ACF);
    
    % Get next frame info
    start=start+frameShift-1;
    stop=start+frameLen-1;
    cnt=cnt+1;
    
    
end
