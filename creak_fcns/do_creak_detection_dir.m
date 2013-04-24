function do_creak_detection_dir(wavDir)

% Function to carry out ANN based creaky voice detection on a directory of
% sound files

%% Initial settings
wavList=dir([wavDir '/*.wav']);
N=length(wavList);

%% Do processing
for n=1:N
   
    basename=regexp(wavList(n).name,'.wav','split');
    basename=char(basename(1));
    disp(basename);
    
    % Load file
    [x,fs]=wavread([wavDir '/' basename '.wav']);
    
    % Do feature extraction and classification
    [prob,dec,t,H2H1,res_p] = CreakyDetection_CompleteDetection2(x,fs);
    
    % Save to .mat
    save([wavDir '/' basename '_creak'],'prob','dec','t','H2H1','res_p')    
    
end