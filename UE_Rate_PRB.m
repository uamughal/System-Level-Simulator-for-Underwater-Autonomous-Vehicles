function [ ue_sinr ue_rate_prb ue_bs_rx_w ue_sinr_eff] = UE_Rate_PRB(apl,ue_bs_rx,ue_bs,N_power,N_prb,N_block,antenna_configuration)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nubsc=sum(apl(:,4)==0);
[nue,nbs]=size(ue_bs);

ue_bs_rx_w=(10.^(ue_bs_rx./10.))/1000.;    %change dBm to W
ue_rx = ue_bs_rx_w.*reshape(repmat(double(ue_bs(:,1:end)),1,N_prb*N_block),nue,nbs,N_prb*N_block);
ue_interf=squeeze(sum(ue_bs_rx_w-ue_rx,2));% nue*nprb
ue_rx=squeeze(sum(ue_rx,2));% nue*nprb

ue_sinr=ue_rx./(ue_interf+N_power);




%update with EESM
ue_rate_prb=(ue_sinr~=0);
[ue_sinr_eff]=EESM_average(ue_sinr,size(ue_sinr,2));
ue_rate_prb=double(ue_rate_prb).*repmat(58.59*TP_LTE_function(ue_sinr_eff,antenna_configuration,'downlink'),1,size(ue_rate_prb,2));
  
end  %--------------------------------end of CL-----------------------------%

