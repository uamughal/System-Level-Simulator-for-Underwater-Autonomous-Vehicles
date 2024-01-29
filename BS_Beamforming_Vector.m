function [w_pha,w_amp]=BS_Beamforming_Vector(M)

n_w=4^(M-1);                    % number of tx_bf code entries
w_pha=zeros(n_w,M);             % phase code book
w_amp=zeros(n_w,M);             % amplitude code book
counter=1;
if M==4
    for n1=1:4
        for n2=1:4
            for n3=1:4
                w_pha(counter,:)=[1 n1,n2,n3];
                counter=counter+1;
            end
        end
    end
elseif M==2
    for n1=1:4
        w_pha(counter,:)=[1 n1];
        counter=counter+1;
    end
elseif M==1
    w_pha=1;
end
% transmit beamforming mode 1
w_pha=exp(1j*pi/2*(w_pha-1))/sqrt(M);
w_amp=ones(n_w,M);