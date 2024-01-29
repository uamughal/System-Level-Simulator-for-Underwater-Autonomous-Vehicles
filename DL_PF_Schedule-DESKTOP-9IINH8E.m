function [ue_tp,ue_scheduled,ue_prb,ue_rate_update,ue_sinr_avg] = DL_PF_Schedule( FD_S,pl_d,ue_layout,bs_layout,ue_bs,ue_rate_prb,ue_rate_avg,TTI_counter,N_prb,N_block,N_prb_mue_exp,ue_bs_rx_w,N_power) 
%UNTITLED2 Summary of this function goes here
% ue_rate ue average rate for 
% Detailed explanation goes here

[nue,nbs]=size(ue_bs);
nmbs=sum(bs_layout(:,4)==0);
ue_scheduled=zeros(1,nue,'int8');
%ue_tp=zeros(nue,N_prb*N_block);
ue_tp=repmat(-inf,nue,N_prb*N_block);
ue_rate_update=ue_rate_avg;
if strcmp(FD_S,'FD_PF')
   
   bs_prb=zeros(nbs,N_block*N_prb);

    for ii=1:nmbs;
        bs_ind=ii;
        ue_ind=find(ue_bs(:,bs_ind)==1); %first column of ue_bs has 12
                                         % non zero element
        if isempty(ue_ind)
            continue
        end
        % ue priority sequence
          % ue_ind=ue_ind(randperm(length(ue_ind)));  %for RR random order
            if sum(ue_rate_avg)==0
              priority_metric=sum(ue_rate_prb(ue_ind,:),2);
           else  priority_metric=sum(ue_rate_prb(ue_ind),2)./ue_rate_avg(ue_ind);
           end
           [priority,ue_order]=sort(priority_metric,'descend');
           ue_ind=ue_ind(ue_order);
           
           prb=repmat(ue_ind.',N_prb_mue_exp,1);
           % ue priority sequence
         
           
          %  prb=[prb,zeros(N_prb_mue_exp,(N_block*N_prb-N_prb_mue_exp*length(ue_ind))/N_prb_mue_exp)];
          % prb=prb(1:N_block*N_prb);% reshape
           prb=reshape(prb,1,(size(prb,1)*size(prb,2)));
           if (length(prb)<N_block*N_prb)
               prb=[prb,zeros(1,N_block*N_prb-length(prb))];
           else prb=prb(1:N_block*N_prb);
           end
           
           prb=prb(randperm(N_block*N_prb));
           ue_scheduled(unique(prb(prb~=0)))=1;
           bs_prb(bs_ind,:)=prb;
           % update the averagerate 
           for jj=1:length(ue_ind)
               ue_id=ue_ind(jj);
               n_prb_ue=sum(prb==ue_id);
               ue_rate_update(ue_id)=((ue_rate_avg(ue_id)*TTI_counter)+ue_rate_prb(ue_id)*n_prb_ue)/(TTI_counter+1);
           end    
      % UE index in each subband and time slot
    %end   
  
    end
      
   end
         

   % number of PRB used by each BS
      n_prb_bs=sum(bs_prb~=0,2);
      n_prb_bs(n_prb_bs==0)=Inf;
     % n_prb_bs(nhbs+1:end)=N_prb;
     % bs_prb(nhbs+1:end,:)=1;
      % ------------ UE received SINR in each PRB ------------------ %
      ue_sinr=zeros(nue,N_prb*N_block);
      % ------------ UE received SINR in each PRB ------------------ %
 %   bs_tx_prb=10.^([repmat(pibs_tx_dBm,1,npibs),repmat(mbs_tx_dBm,1,nmbs)]/10)/1000./n_prb_bs.';
  %if ABS_indicator==0
  % bs_tx_prb=10.^([repmat(pibs_tx_dBm,1,npibs),repmat(mbs_tx_dBm,1,nmbs)]/10)/1000./N_prb;
  %else
  %  bs_tx_prb=[10.^(repmat(pibs_tx_dBm,1,npibs)/10)/1000./N_prb,repmat(0,1,nmbs)]; 
  %end
    

  %  rx=reshape(repmat(bs_tx_prb,nue,N_block*N_prb),nue,nbs,N_prb*N_block)./pl;
    rx=ue_bs_rx_w;% considering ABS indicator
   % tx_ind=bs_prb~=0;
   % rx=rx.*shiftdim(reshape(repmat(tx_ind,1,nue),nbs,N_prb*N_block,nue),2); %only eNBs with users will Tx

    serve_ind=reshape(repmat(ue_bs,1,N_block*N_prb),nue,nbs,N_prb*N_block);
    sig=squeeze(sum(rx.*double(serve_ind),2));
    int=squeeze(sum(rx,2))-sig;
    sinr=sig./(int+N_power);
      % SINR - throughput mapping
    active_ue=find(ue_scheduled~=0);
    ue_prb=zeros(nue,N_prb*N_block);
    for ii=1:length(active_ue)
    ue=active_ue(ii);
    [bs_ind prb_ind]=find(bs_prb==ue);
%     sinr_ind=sub2ind([nue N_prb*N_block],repmat(ue,length(prb_ind),1),prb_ind);
    %         ue_tp(ii)=sum(TP_LTE_function(sinr(sinr_ind),2,'downlink')*180);
%     ue_tp(ii,prb_ind)=TP_LTE_function(sinr(sinr_ind),2,'downlink')*180;
    ue_tp(ue,prb_ind)=10*log10(sinr(ue,prb_ind));
    ue_prb(ue,prb_ind)=1;
    ue_sinr(ue,prb_ind)=sinr(ue,prb_ind);
    end
    ue_sinr_sum=sum(ue_sinr,2);
    ue_prb_sum=sum(ue_prb,2);%total prb for each UE
    ind=ue_sinr_sum==0;
    ue_sinr_sum(ind)=[];
    ue_prb_sum(ind)=[];
    ue_sinr_avg=ue_sinr_sum./ue_prb_sum;


 
end
 
