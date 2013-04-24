function [NAQ,QOQ,H1H2,HRF] = get_NAQ_QOQ_H1H2(glot,fs,GCI)

% Function to calculate NAQ, QOQ and H1-H2 values from an estimated 
% glottal source signal. All parameter values are aligned to the inputted
% GCIs.

%% Initial settings
F0min=20;
F0max=500;
NAQ=zeros(1,length(GCI));
QOQ=zeros(1,length(GCI));
H1H2=zeros(1,length(GCI));
HRF=zeros(1,length(GCI));

qoq_level=0.5; 
T0_num=3;

glot_int=integrat(glot,fs);
glot_shift=0.5/1000*fs;

%% Do processing
for n=1:length(GCI)
    
    % Get glottal pulse compensated for zero-line drift
    if n==1
        start=1;
        stop=GCI(n);
        T0=GCI(n+1)-GCI(n);
    else start=GCI(n-1);
        stop=GCI(n);
        T0=GCI(n)-GCI(n-1);
    end
    F0=fs/T0;

    if isinf(F0)==0 && T0~=0 && F0 > F0min && F0<F0max 
    
        glot_int_comb=[glot_int(start) glot_int(stop)];
        if start~=stop
            if length(glot_int_comb)>1
                line=interp1(linspace(1,stop-start+1,2),glot_int_comb, ...
                    1:stop-start+1);
            else line=glot_int_comb;
            end
        else line=0;
        end
        glot_int_cur=glot_int(start:stop);
        glot_int_comp=glot_int_cur-line;

        if stop+glot_shift <= length(glot)
            stop2=stop+glot_shift;
        else stop2=stop;
        end
        glot_cur=glot(start:stop2);

        % Get frame positions for H1-H2 parameter
        if GCI(n)-round((T0*T0_num)/2) > 0
            f_start = GCI(n)-round((T0*T0_num)/2);
        else f_start=1;
        end
        if GCI(n)+round((T0*T0_num)/2) <= length(glot)
            f_stop = GCI(n)+round((T0*T0_num)/2);
        else f_stop=length(glot);
        end
        f_frame=glot(f_start:f_stop);
        f_win=f_frame(:).*hamming(length(f_frame));
        f_spec = 20*log10(abs(fft(f_win,fs)));

        % Get H1-H2 measurements
        [h_idx,h_amp]=findpeaks(f_spec,[],F0/2);
        if length(h_idx) < 5
            H1H2(n)=0;
        else 
            [~,f0_idx]=min(abs(h_idx-F0)); % Find closest peak to F0
            [~,f02_idx]=min(abs(h_idx-(F0*2))); % Find closest peak to F0*2
            H1H2(n)=h_amp(f0_idx)-h_amp(f02_idx);
            HRF(n) = sum(h_amp(f0_idx))/h_amp(f0_idx(1));
        end


        % Get NAQ, QOQ and H1-H2 measurements
        d_peak = max(abs(glot_cur));
        [f_ac,max_idx] = max(glot_int_comp);
        Amid=f_ac*qoq_level;
        [T1,T2] = findAmid_t(glot_int_comp,Amid,max_idx);

        NAQ(n) = (f_ac/d_peak)*F0;
        QOQ(n)=(T2-T1)/(fs/F0);

        % Plots
    %     subplot(211)
    %     glot_int_compNorm=glot_int_comp/max(glot_int_comp);
    %     plot(glot_int_compNorm), hold on, 
    %     plot(glot_cur/max(abs(glot_cur)),'r'), 
    %     plot(T1,glot_int_compNorm(T1),'rx'),plot(T2,glot_int_compNorm(T2),'rx')
    %     hold off, 
    %     subplot(212)
    %     plot(f_spec),hold on, plot(h_idx,h_amp,'rx'), hold off
    %     xlim([0 fs/2])
    %     [f0_idx f02_idx]
    %     pause(0.5);
    end
end

QOQ(QOQ<0|QOQ>1)=0;
        

