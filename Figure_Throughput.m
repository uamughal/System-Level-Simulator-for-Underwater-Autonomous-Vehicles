% function [hue5 hue50 mue_in5 mue_in50 mue_out5 mue_out50]=Figure_Throughput(UE_TP,UE_schedule,UE_type,UE_PRB,TP_type)
function [mue_out_sinr cue_out_sinr reue_out_sinr ue_out_sinr ue_out_efficiency mue_pos cue_pos reue_pos]=Figure_Throughput(UE_TP,UE_schedule,UE_type,UE_PRB,TP_type,Antenna_configuration,ue_type,ue_layout)

hue_tp=[];
mue_in_tp=[];
mue_out_tp=[];
ue_sinr=[];
n=size(UE_TP,1);
ue_out_efficiency=[];
mue_sinr=[];
%PUE
cue_sinr=[];
reue_sinr=[];
%UE_pos
mue_pos=[];
cue_pos=[];
reue_pos=[];

for ii=1:n
    sinr_dB=UE_TP{ii};
  %  sinr_dB(sinr_dB==0)=-Inf;
  
   % tp=TP_LTE_function(10.^(sinr_dB/10),[],'uplink')*180;

   tp=TP_LTE_function(10.^(sinr_dB/10),Antenna_configuration,'downlink')*58.59;
   ind_tp=abs(tp)==Inf;
   tp(ind_tp)=0;
    schedule=UE_schedule{ii};
    type=UE_type{ii};
    prb=UE_PRB{ii};
    
    % user throughput
    if strcmp(TP_type,'user')
        hue_tp=[hue_tp; sum(tp(schedule==1 & type==1,:),2)];
        mue_in_tp=[mue_in_tp; sum(tp(schedule==1 & type==2,:),2)];
        mue_out_tp=[mue_out_tp; sum(tp(schedule==1 & type==3,:),2)];
        
        ue_sinr=[ue_sinr; sinr_dB(schedule==1,:)]; % To find the overall SINR of the users
        ue_out_sinr= max(ue_sinr,[],2); %why max//
        ue_out_efficiency=[ue_out_efficiency;TP_LTE_function(10.^(ue_out_sinr/10),Antenna_configuration,'uplink')];
        
        mue_sinr=[mue_sinr;sinr_dB(schedule==1&ue_type==1,:)];
        mue_out_sinr= max(mue_sinr,[],2); %why max//
        
        cue_sinr=[cue_sinr;sinr_dB(schedule==1&ue_type==2,:)];
        cue_out_sinr= max(cue_sinr,[],2); %why max//
        
        reue_sinr=[reue_sinr;sinr_dB(schedule==1&ue_type==3,:)];
        reue_out_sinr= max(reue_sinr,[],2); %why max//
        
        %position
        mue_pos=[mue_pos;ue_layout(schedule==1&ue_type==1,1:2)];
        cue_pos=[cue_pos;ue_layout(schedule==1&ue_type==2,1:2)];
        reue_pos=[reue_pos;ue_layout(schedule==1&ue_type==3,1:2)];
        
        
        
    elseif strcmp(TP_type,'PRB')
        hue=tp(schedule==1 & type==1,:);
        mue_in=tp(schedule==1 & type==2,:);
        mue_out=tp(schedule==1 & type==3,:);
        hue_prb=prb(schedule==1 & type==1,:);
        mue_in_prb=prb(schedule==1 & type==2,:);
        mue_out_prb=prb(schedule==1 & type==3,:);
        
        hue_tp=[hue_tp; hue(hue_prb==1)];
        mue_in=mue_in(mue_in_prb==1);
        mue_in_tp=[mue_in_tp; mue_in(:)];
        mue_out_tp=[mue_out_tp; mue_out(mue_out_prb==1)];
    end
end


%% cdf figure f(x)

%if ~isempty(hue_tp)
%   [f_hue,x_hue]=ecdf(hue_tp);
%   [val hue5_ind]=min(abs(f_hue-0.1));
%   hue5=x_hue(hue5_ind);
%   [val hue50_ind]=min(abs(f_hue-0.5));
%   hue50=x_hue(hue50_ind);
%else
%    x_hue=[];
%    f_hue=[];
%    hue5=0;
%    hue50=0;
%end
%if ~isempty(mue_in_tp) 
%    [f_mue_in,x_mue_in]=ecdf(mue_in_tp);
%    [val mue_in5_ind]=min(abs(f_mue_in-0.1));
%    mue_in5=x_mue_in(mue_in5_ind);
%    [val mue_in50_ind]=min(abs(f_mue_in-0.5));
%    mue_in50=x_mue_in(mue_in50_ind);
%else
%    x_mue_in=[];
%    f_mue_in=[];
%    mue_in5=0;
%    mue_in50=0;
%end
%if ~isempty(mue_out_tp)
%    [f_mue_out,x_mue_out]=ecdf(mue_out_tp);
%    [val mue_out5_ind]=min(abs(f_mue_out-0.1));
%    mue_out5=x_mue_out(mue_out5_ind);
%    [val mue_out50_ind]=min(abs(f_mue_out-0.5));
%   mue_out50=x_mue_out(mue_out50_ind);
%    [f_ue_out_sinr,x_ue_out_sinr]=ecdf(ue_out_sinr);
    
%else
%    x_mue_out=[];
%    f_mue_out=[];
%    mue_out5=0;
%    mue_out50=0;
%end