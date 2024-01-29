function [ue_bs_rx] = UE_BS_Rx(ue_layout,bs_layout,mbs_tx_dBm,pl_d,pl_slow,nprb_m,N_block)
%UNTITLED Summary of this function goes here
%this can only be used when there is no henb&hue
%   Detailed explanation goes here

nbs=size(bs_layout,1);
nmbs=nbs;
n_ue=size(ue_layout,1);

%ue_bs_rx_bias=zeros(n_ue,nmbs+npibs);
%ue_bs_rx=zeros(n_ue,nmbs+npibs);
mbs_tx_dBm_prb=mbs_tx_dBm-10*log10(nprb_m);
tx_mbs=repmat(mbs_tx_dBm_prb,n_ue,nmbs);
% tx_ue_bs_bias=[tx_mbs];
tx_ue_bs=[tx_mbs];

% ue_bs_rx_bias=reshape(repmat(tx_ue_bs,1,N_block*nprb_m),n_ue,nmbs,N_block*nprb_m)-reshape(repmat(pl_slow,1,nprb_m*N_block),n_ue,nbs,nprb_m,N_block);
% ue_bs_rx_bias=sum((10.^(ue_bs_rx_bias./10.)),3)./(N_block*nprb_m);%mw
% ue_bs_rx_bias=10*log10(ue_bs_rx_bias);%dBm
%ue_bs_rx_avg_dBm=10.*log10(ue_bs_rx_bias)-[repmat(10,n_ue,npibs),zeros(n_ue,nmbs)];
ue_bs_rx=reshape(repmat(tx_ue_bs,1,N_block*nprb_m),n_ue,nmbs,N_block*nprb_m)-pl_d(:,1:end,:);
end

