function plot_spectogram(x,fs,window,nfft)

if nargin < 3
    window=256;
end
if nargin < 4
    nfft=1024;
end

% calculate the table of amplitudes
[B,f,t] = specgram(x,nfft,fs,window,round(window*.8));
%[B,f,t]=specgram(x,1024,fs,256,192);

%
% calculate amplitude 50dB down from maximum
bmin=max(max(abs(B)))/8000;
%
% plot top 50dB as image
imagesc(t,f,20*log10(max(abs(B),bmin)/bmin));
%
% label plot
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%
% build and use a grey scale
lgrays=zeros(100,3);
for i=1:100
    lgrays(i,:) = 1-i/100;
end
colormap(lgrays);