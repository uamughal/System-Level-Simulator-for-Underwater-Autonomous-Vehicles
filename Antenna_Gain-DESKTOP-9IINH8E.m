function ant_gain = Antenna_Gain( MBS_ant,MBS_div_gain,w_pha,w_amp,UBSC_layout,UBS_layout,ue_bs,UBSC_ff,N_block,N_prb )

nbs=size(UBSC_layout,1);
nue=size(UBS_layout,1);
n_ant=size(w_pha,2);
hue_ind=UBS_layout(:,7)~=0;
if ~isempty(UBSC_ff)
    n_block=size(UBSC_ff,4);
    n_prb=size(UBSC_ff,3);
else
    n_block=1;
    n_prb=1;
end
% add
   n_block=N_block;
    n_prb=abs(N_prb);

% find UE-MBS link DOA
comp_sep=repmat(UBS_layout(:,1)+1j*UBS_layout(:,2),1,nbs)-repmat((UBSC_layout(1:end,1)+1j*UBSC_layout(1:end,2)).',nue,1);
ue_doa=angle(comp_sep)-angle(repmat(UBSC_layout(1:end,4).',nue,1));   %angle(H) returns the phase angles, in radians, of a matrix with
                                                                      %complex elements.
[index1,index2]=find(ue_doa>pi);
ue_doa=ue_doa+sparse(index1,index2,-2*pi,nue,nbs);
[index1,index2]=find(ue_doa<-pi);
ue_doa=ue_doa+sparse(index1,index2,2*pi,nue,nbs);
ue_doa=180/pi*ue_doa; 

% calculate UE-MBS antenna gain
ant_gain=zeros(nue,nbs,n_prb,n_block);   %gain is calculated for every PRB
if strcmp(MBS_ant,'Omi')
%     ant_gain(:,nhbs+npibs+1:nbs,:,:)=reshape(repmat(14-min(12*(ue_doa(:)/70).^2,20),n_block*n_prb,1),nue,nbs-nhbs-npibs,n_prb,n_block)+MBS_div_gain;
    ant_gain(:,1:nbs,:,:)=reshape(repmat(5,n_block*n_prb*nue*(nbs),1),nue,nbs,n_prb,n_block);
%     ant_gain(:,nhbs+1:nhbs+npibs,:,:)=reshape(repmat(5,n_block*n_prb*nue*npibs,1),nue,npibs,n_prb,n_block);
end
end