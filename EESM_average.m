function [ ue_bs_sinr_eff ] = EESM_average( ue_prb_sinr,N_prb )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% unit for ue_bs_prb_sinr 1
%Modulation	Code Rate ¦Â
% QPSK	      1/3	 1.49
% 	          1/2	 1.57
% 	          2/3  	 1.69
% 	          3/4	 1.69
% 	          4/5	 1.65
% 
% 
% 16QAM       1/3	 3.36
% 	          1/2	 4.56
% 	          2/3	 6.42
% 	          3/4	 7.33
% 	          4/5	 7.68

%------------Link Adaptation--------------------------%
    
%% MCS selection 
%------------------MCS_beta----------------%
nue=21;
      cqi_index=zeros(nue,1);
  for ue=1:nue 
          
       if (ue_prb_sinr(ue)<-8.3)
          cqi_index(ue)=1;
      elseif (ue_prb_sinr(ue)>=-6.2)&&(ue_prb_sinr(ue)<-6.1)
          cqi_index(ue)=2;
      elseif (ue_prb_sinr(ue)>=-6.1)&&(ue_prb_sinr(ue)<-5.6)
          cqi_index(ue)=3;  
      elseif (ue_prb_sinr(ue)>=-5.6)&&(ue_prb_sinr(ue)<-4.5)
          cqi_index(ue)=4;  
      elseif (ue_prb_sinr(ue)>=-4.5)&&(ue_prb_sinr(ue)<-3.6)
          cqi_index(ue)=5;  
      elseif (ue_prb_sinr(ue)>=-3.6)&&(ue_prb_sinr(ue)<-2.9)
          cqi_index(ue)=6;  
      elseif (ue_prb_sinr(ue)>=-2.9)&&(ue_prb_sinr(ue)<-6.1)
          cqi_index(ue)=7;  
      elseif (ue_prb_sinr(ue)>=-6.1)&&(ue_prb_sinr(ue)<0.6)
          cqi_index(ue)=8;  
      elseif (ue_prb_sinr(ue)>=0.6)&&(ue_prb_sinr(ue)<1.9)
          cqi_index(ue)=9;  
      elseif (ue_prb_sinr(ue)>=1.9)&&(ue_prb_sinr(ue)<2.9)
          cqi_index(ue)=10;  
      elseif (ue_prb_sinr(ue)>=2.9)&&(ue_prb_sinr(ue)<3.2)
          cqi_index(ue)=11; 
      elseif (ue_prb_sinr(ue)>=3.2)&&(ue_prb_sinr(ue)<4.1)
          cqi_index(ue)=12;  
      elseif (ue_prb_sinr(ue)>=4.1)&&(ue_prb_sinr(ue)<5.8)
          cqi_index(ue)=13;  
      elseif (ue_prb_sinr(ue)>=5.8)&&(ue_prb_sinr(ue)<10.4)
          cqi_index(ue)=14;  
      elseif ue_prb_sinr(ue)>=10.4
        cqi_index(ue)=15;
      end
      cqi_index(cqi_index==0)=1;
  [beta_LA]=find_beta_MCSselection(cqi_index,nue);  %find the beta first(choosing MCS according to cqi index)
  ue_bs_sinr_eff=-beta_LA(ue)*log(sum(exp(-ue_prb_sinr./beta_LA(ue)),2)./N_prb);  %actual EESM formula
  end 

% beta_table=[2.1,4.5,6.7];
% beta=beta_table(1);
% ue_bs_sinr_eff=-beta*log(sum(exp(-ue_prb_sinr./beta),2)./N_prb);

end

