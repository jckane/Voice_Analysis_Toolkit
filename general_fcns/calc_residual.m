function [vector_res,ar_lpc,e_lpc]=calc_residual(x,x_lpc,ord_lpc,GCI)

% Function to carry out LPC analysis and inverse filtering, used in IAIF.m 
% function.

% USAGE:    
%       Input:
%             x       : signal to be inverse filtered
%             x_lpc   : signal to carry out LPC analysis on
%             ord_lpc : LPC prediction order
%             GCI     : Glottal closure instants (in samples)
%
%       Output:
%             vector_res : residual after inverse filtering
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Coded by John Kane @ The Phonetics and Speech Lab %%%%%%%%%%%
%% Trinity College Dublin, August 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial settings
vector_res=zeros(1,length(x)); % allocate space
ze_lpc=zeros(1,ord_lpc); % "

ar_lpc=zeros(ord_lpc+1,length(GCI));
e_lpc=zeros(ord_lpc+1,length(GCI));
%plot(x_lpc)

%% Do processing
for n=1:length(GCI)
    
    % Get framing information
    if n > 1
        T0_cur=GCI(n)-GCI(n-1);
    else T0_cur=GCI(n+1)-GCI(n);
    end

    if GCI(n)-T0_cur > 0 && GCI(n)+T0_cur <= length(x)
        
        start=GCI(n)-T0_cur;
        stop=GCI(n)+T0_cur;
        
        % Do LPC analysis
        frame_lpc=x_lpc(start:stop);
        if length(frame_lpc) > ord_lpc*1.5
            frame_wind=frame_lpc.*hamming(length(frame_lpc))';    
            [ar,e]=lpcauto(frame_wind,ord_lpc);     
            ar=real(ar);

            ar_lpc(:,n)=ar(:);
            e_lpc(n)=e;


            % Do inverse filtering
            if n > 1 && exist('frame_res','var')

                last_input=fliplr(frame_res);
                last_output=fliplr(residual);
                ze_lpc=filtic(ar,sqrt(e),last_output,last_input);
            end

            frame_res=x(start:stop);

            % calculation of the LPC residual
            try
                [residual]=filter(ar,sqrt(e),frame_res,ze_lpc);
            catch residual =frame_res;
                disp('Problem with IAIF filtering')
            end
            residual_win=residual(:).*hamming(length(residual));

          %  plot(frame_wind), hold on,plot(residual_win,'r'), hold off
          %  pause(.3)

            vector_res(start:stop)=vector_res(start:stop)+residual_win';
        end
    end

       
end