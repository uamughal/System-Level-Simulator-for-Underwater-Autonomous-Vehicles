%%%------------------- OSTN Downlink SLS Link 1-------------------------%%

clear all;
close all;

ISD = 40000;        % inter site distance 40km
n_rings=1;          % 1 tier
sim_resolution=10 ;  % how far it has to jump to cover the grid area
UBS_free=10;        % this parameter is not in use
D1=1300;            % Distances of UBS's,
D2=5800;            % We should use these values
D3=11500;           % Otherwise it doesn't give correct values
FD_S='FD_PF';       % Scheduling(use the same)
N_prb_mue_exp=10;   % expected number of PRB per UE
load FixedData_SLS;
CL=1;               % closed loop(feedback issues)


%--------------- Frame Structure(no change)-----------------------------------
 fc = 5.5e3;
 Num_OFDMsym_per_Frame = 7*8;
     R = 5e3;              % filter bw
        Len_PRB_freq = 62; % Number of subcarriers per PRB (496/8)
        Nfft = 512;        % fft size
        Ncp = 113;         % cp length=113 samples
        gb1 = 6;           % guard band(left)
        gb2 = 7;           % guard band(right)
     
        fs_baseband = 5e3; % baseband freq.
Num_Pilot_Spacing = 1;
codingrate = 1/2;
  Len_PRB_time = 6;
  Num_PreambleSym = 2;
Num_DataOFDMsym_per_Frame = Num_OFDMsym_per_Frame - Num_PreambleSym;
Num_DC_subcarrier = 3;
Nefft = Nfft - gb1 - gb2 - Num_DC_subcarrier; 
pilot_space = 2;
Num_PilotSC_per_PRB = length(1:pilot_space:Len_PRB_freq);
Num_DataSC_per_PRB = Len_PRB_freq - Num_PilotSC_per_PRB;
idx_Pilot_PRB = (1:pilot_space:Len_PRB_freq).';              %.(dot) means all, ' means transpose
idx_DataSC_PRB = setdiff((1:Len_PRB_freq), idx_Pilot_PRB).'; % C = setdiff(A,B) for vectors A and B, returns the values in A that 
                                                             % are not in B with no repetitions. C will be sorted.
idx_dc = round(Nfft/2)-1:round(Nfft/2)+1;
idx_UsefulSC_IFFT = [gb1+1:idx_dc(1,1)-1 idx_dc(1,3)+1:Nfft-gb2];

% switch Modorder
%     case 2
%         MODULATION = 'bpsk';
%     case 4
%         MODULATION = 'qpsk';
%     case 16
%         MODULATION = '16qam';
% end

CRC_Len_Bits = 24;                % bits
Len_PRB_freq_per_CW = 1;
Len_PRB_time_per_CW = 1;
Num_RequiredPRB_per_CW = 1;


% Num_Sym_Per_CW = Len_PRB_freq_per_CW*Num_DataSC_per_PRB*Len_PRB_time_per_CW*Len_PRB_time;
% Num_Bits_Per_CW = Num_Sym_Per_CW*log2(Modorder)*codingrate - CRC_Len_Bits;
% msg_Len = Num_Bits_Per_CW;

% TxPilot = TxPilot(1:length(idx_Pilot_IFFT),1);

xx = xx_st{1,11};
TxFrame = xx;
% yy = reshape( TxFrame, Nfft+Ncp, Num_DataOFDMsym_per_Frame );
% yy_noCP = yy(Ncp+(1:Nfft), :);

%%---------------upto here is frame structure-----------------%%


%%%--------------------UBSC Parameters-------------------------%%
%N_prb=TxFrame;

N_prb=83;                            %number of physical resource block
N_block=1;                           %number of blocks
UBSC_tx_dBm=50;                      %transmit power in dbm
antenna_configuration=1;             %SISO
UBSC_ant_type='omi';
M=1;                                 % number of antennas
MBS_ant_type='Omi'; 
MBS_div_gain=0;

%%-----------------Environment Parameters---------------------%%

Wavelength=(5*10^3)/(5.5*10^3);            % wavelength
BW_PRB=58.59;                              % PRB spectrum width
fkz=5.5*10^3;                              % Carrier Frequency
BW=5*10^3;                                 % Bandwidth
 

% %% Path Loss
% xa = fkz.^2;
% fa = (((0.11*(xa))./(1+(xa))) + ((44*(xa))./(4100+(xa))) + (2.75*(xa)/10000) + 0.003);
% k=1.5;     % Spreading Factor
%%pl=((k*log10(d))+(d*log10(fa)));
%%pl_a=10.^(pl/10);

%%---------------------------Ambient Noise--------------------------%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Turbulance Noise%%%%%%%%%%%%%%%%%%%%%%%%%%

N_t=17-30*log10(fkz);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Shipping Noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=0;                  % shipping activity factor (value 0~1 =low~high)
N_s=40+(20*(s-0.5))+26*log10(fkz)-60*log10(fkz+0.03);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Wind Noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=0;              % wind speed in m/s
N_w=50+(7.5*(w^0.5))+20*log10(fkz)-40*log10(fkz+0.4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Thermal Noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_th=-15+20*log10(fkz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Total Noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_power=10*log10(10.^(N_t/10)+10.^(N_s/10)+10.^(N_w/10)+10.^(N_th/10));    %dB Value
N_power_a=(10.^(N_t/10)+10.^(N_s/10)+10.^(N_w/10)+10.^(N_th/10));                    %Watt Value


N_prb_UBS_exp=10;                % expected number of PRB per UBS
is_fast_fading=1;

%%---------------------- Simulation Parameters--------------------%%
isplot=1;
drops=1;               %number of drops means how many times we calculate the value
TTI_MAX=50;            % TTI= round trip between UE and BS, 
                       % so 50 means we will calculate for 50 round trips.
                       % Thus it will give us some data to calculate the average.
                       % if too less "tti" then it's not enough data to calculate the average.
Frame_structure='TDD_1';

%%-----------------------UBSC Layout--------------------------%%
[UBSC_layout UBS_grid grid_cell_ind]=UBnew_Layout(ISD,n_rings,isplot,sim_resolution,UBS_free);
n_cells=size(UBSC_layout,1);     % total no. of sectors
N_UBS=3;                         % number of UBS

%%%%%%%%%%%%%%%%%%%BS transmit antenna weight vector%%%%%%%%%%%%%%%%%%%%
[w_pha,w_amp]=BS_Beamforming_Vector(M); % w_pha=phase code book,
                                        % w_amp=amplitude code book.
                                        %antenna weight vector= it makes
                                        %the omnidirectional antenna to a
                                        %directional antenna.
ue_rate_avg=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%initilize result collector%%%%%%%%%%%%%%%%%%%%%

UE_TP=cell(drops,1);   %cell arrays contain either list of text, combination
                       %of text and numbers, or numeric arrays of different
                       %sizes
UE_type=cell(drops,1); 
UE_scheduled=cell(drops,1);
UE_BS=cell(drops,1);
PRB=cell(drops,1);
rate_avg_collector=[];
rmm_e=[];
r11_e=[];
r12_e=[]; 
UE_Sinr=[];
UE_sinr_collector=[];
UE_thp_collector=[];
UE_efficiency_collector=[];
MUE_pos_collector=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Main Loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:drops
    ii
    area_grid=UBS_grid;
    cell_ind=grid_cell_ind;

%UBS Layout
%[UBS_layout]=UBS_Layout(UBSC_layout,isplot,0,'fixed',area_grid,cell_ind,N_UBS,n_cells,ISD,D1,D2,D3);
    [UBS_layout]=UBS_Layout(UBSC_layout,D1,D2,D3);
    
  %%%%%%%%%%%%%%%%%%%results collector for each drop%%%%%%%%%%%%%%%%%%%%%

UE_TP=cell(TTI_MAX,1);
UE_type=cell(TTI_MAX,1); 
UE_scheduled=cell(TTI_MAX,1);
UE_BS=cell(TTI_MAX,1);
PRB=cell(TTI_MAX,1);

[UBS_layout UBSC_layout UBS_UBSC]=BS_Association(UBS_layout,UBSC_layout);
ue_rate_avg=zeros(size(UBS_layout,1),1);
antenna_gain=Antenna_Gain(MBS_ant_type,MBS_div_gain,w_pha,w_amp,UBSC_layout,UBS_layout,UBS_UBSC,[],N_block,N_prb);
 [DL_UL_indicators]=TTI_indicator(Frame_structure);  % frame configuration
 
  for TTI_counter=1:TTI_MAX
        TTI_counter

%%%%%%%%%%%%%%%%%%%%%%%%%%%subframe identification%%%%%%%%%%%%%%%%%%%%%
    
subframe_counter=mod(TTI_counter-1,10)+1;     %Modulus-> takes positive
                                              %values. in order to avoid
                                              %negative values we used mod.
                                              %another term for absolute
                                              %value.
 
%%%%%%%%%%%%%%%%%%%%%%Not working%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DL_UL_indicator=DL_UL_indicators(subframe_counter);
    
%%-----------------------------Path Loss-------------------------------%%    
    pl_avg=Path_Loss(UBSC_layout,UBS_layout);
   %pl_avg=reshape(pl_avg1,21,7);
   

%%--------------OSTN Channel Model------------------------%% 
    pl=fd_ch(UBSC_layout,UBS_layout);   %fading channel model

    mbs_ffa=10*log10(abs(pl));
    mbs_ffa(isinf(mbs_ffa))=0; %For example, isinf([pi NaN Inf -Inf])--->
                               %is [0 0 1 1]. 1 for infinite number, 0 for
                               %others.
    nue=size(UBS_layout,1);    %UE stuffs
    nbs=size(UBSC_layout,1);   %UBSC stuffs     
    mbs_ff=reshape(repmat(mbs_ffa,1,N_prb*N_block),nue,nbs,N_prb*N_block);
    pl_slow=pl_avg;                              
    pl_slow(:,1:end)=pl_slow(:,1:end)-antenna_gain(:,1:end,1);        
% pl_slow(:,nhbs+1:end)=pl_slow(:,nhbs+1:end);                 
    pl_d=reshape(repmat(pl_slow,1,N_prb*N_block),nue,nbs,N_prb,N_block);              
    pl_d(:,1:end,:)=pl_d(:,1:end,:)-mbs_ff;                                 
    pl_d(:,1,:)=pl_d(:,1,:)-antenna_gain(:,1,:);
  
   % determine MUE serving BS by selecting BS with least path loss
   [UBS_UBSC_Rx]=UE_BS_Rx(UBS_layout,UBSC_layout,UBSC_tx_dBm,pl_d,pl_slow,N_prb,N_block);
   % RSRP with bias for MBS&PicoBS dimension n_ue*(n_mbs+n_pibs)
    
   if TTI_counter==1
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% UBS and UE association%%%%%%%%%%%%%%%%%%%%%
   
UBS_ue=find(UBS_layout(:,5)==0);
   
[val,serving_bs]=max(UBS_UBSC_Rx(UBS_ue,:),[],2);  %maximum component
                       %Search google keep for----->>Sim_Solution
                       %val---> maximum comp., serving_bs---> row no.
   
   mpbs_ind=sub2ind(size(UBS_UBSC),UBS_ue,serving_bs);   %IND = sub2ind(SIZ,I,J) returns the linear index equivalent to the
                                                         %row and column subscripts in the arrays I and J for a matrix of
                                                         %size SIZ.
                                                         %linear index=postion of
                                                         %the value in actual matlab
                                                         %storage.
                                                         
                                                         %subscript=indexing of row
                                                         %and coloumn. 
                                                         
   
   UBS_UBSC(mpbs_ind)=1;
%[ ue_type ] = UE_Classification(UBS_UBSC,UBSC_layout,ue_bs_Rx,N_prb,N_block,N_power,SINR_threshold);
   end
[ue_sinr,ue_rate_prb,ue_bs_rx_w,ue_sinr_eff]=UE_Rate_PRB(UBSC_layout,UBS_UBSC_Rx,UBS_UBSC,N_power,N_prb,N_block,antenna_configuration);

%------------------------- DL scheduling -----------------------%
if DL_UL_indicator==1
if strcmp(FD_S,'FD_PF') % TF = strcmp(S1,S2) compares S1 and S2 and 
                        % returns logical 1 (true)
                        % if they are identical, and returns
                        %logical 0 (false) otherwise.
   [ue_tp,ue_scheduled,prb_allocation,ue_rate_update,ue_sinr_avg]=DL_PF_Schedule(FD_S,pl_d,UBS_layout,UBSC_layout,UBS_UBSC,ue_rate_prb,ue_rate_avg,TTI_counter,N_prb,N_block,N_prb_mue_exp,ue_bs_rx_w,N_power);
else if strcmp(FD_S,'FD_RR')  
   [ue_tp,ue_scheduled,prb_allocation,ue_rate_update,ue_sinr_avg]=DL_RR_Schedule(FD_S,pl_d,UBS_layout,UBSC_layout,UBS_UBSC,ue_rate_prb,ue_rate_avg,TTI_counter,N_prb,N_block,N_prb_mue_exp,ue_bs_rx_w,N_power);    
    end
end
  ue_rate_avg=ue_rate_update;
else ue_rate_avg=ue_rate_avg*TTI_counter/(TTI_counter+1);
end




%%-------------------------Link Adaptation--------------------------%
% 
% if CL==1 
%     
%% MCS selection     
%---------------------------MCS_beta_offset-------------------------%

      cqi_index=zeros(nue,1);
  for ue=1:nue 
   
      if (ue_sinr_eff(ue)<-8.3)
          cqi_index(ue)=1;
      elseif (ue_sinr_eff(ue)>=-6.2)&&(ue_sinr_eff(ue)<-6.1)
          cqi_index(ue)=2;
      elseif (ue_sinr_eff(ue)>=-6.1)&&(ue_sinr_eff(ue)<-5.6)
          cqi_index(ue)=3;  
      elseif (ue_sinr_eff(ue)>=-5.6)&&(ue_sinr_eff(ue)<-4.5)
          cqi_index(ue)=4;  
      elseif (ue_sinr_eff(ue)>=-4.5)&&(ue_sinr_eff(ue)<-3.6)
          cqi_index(ue)=5;  
      elseif (ue_sinr_eff(ue)>=-3.6)&&(ue_sinr_eff(ue)<-2.9)
          cqi_index(ue)=6;  
      elseif (ue_sinr_eff(ue)>=-2.9)&&(ue_sinr_eff(ue)<-6.1)
          cqi_index(ue)=7;  
      elseif (ue_sinr_eff(ue)>=-6.1)&&(ue_sinr_eff(ue)<0.6)
          cqi_index(ue)=8;  
      elseif (ue_sinr_eff(ue)>=0.6)&&(ue_sinr_eff(ue)<1.9)
          cqi_index(ue)=9;  
      elseif (ue_sinr_eff(ue)>=1.9)&&(ue_sinr_eff(ue)<2.9)
          cqi_index(ue)=10;  
      elseif (ue_sinr_eff(ue)>=2.9)&&(ue_sinr_eff(ue)<3.2)
          cqi_index(ue)=11; 
      elseif (ue_sinr_eff(ue)>=3.2)&&(ue_sinr_eff(ue)<4.1)
          cqi_index(ue)=12;  
      elseif (ue_sinr_eff(ue)>=4.1)&&(ue_sinr_eff(ue)<5.8)
          cqi_index(ue)=13;  
      elseif (ue_sinr_eff(ue)>=5.8)&&(ue_sinr_eff(ue)<10.4)
          cqi_index(ue)=14;  
      elseif ue_sinr_eff(ue)>=10.4
        cqi_index(ue)=3;
      end
    cqi_index(cqi_index==0)=1;
  [beta_LA]=find_beta_MCSselection(cqi_index,nue);
  end 
%   
% end  %--------------------------------end of CL-----------------------------%
 
 
   
 
end 

UE_Sinr=[];
UE_sinr_collector=[];
MUE_sinr_collector=[];
no_mue=0;
MUE_pos_collector=[];

UE_Sinr=[UE_Sinr;ue_sinr_avg];

   
%-------------------------- throughput collect ----------------------- %
   temp=zeros(1,0,'int8');
   temp( UBS_layout(:,8)==0)=3; % 
   temp( UBS_layout(:,8)~=0)=2; % 
   temp(UBS_layout(:,8)<=0 & UBS_layout(:,7)~=0)=1;
    
    ind_remove=find(temp==0);
%     ue_tp(ind_remove)=[];
    ue_tp(ind_remove,:)=[];
    ue_scheduled(ind_remove)=[];
    temp(ind_remove)=[];
    UBS_layout(ind_remove,:)=[];
    prb_allocation(ind_remove,:)=[];
    
    UE_TP{TTI_counter}=ue_tp;
    UE_scheduled{TTI_counter}=ue_scheduled;
    UE_type{TTI_counter}=temp;
    UE_BS{TTI_counter}=UBS_UBSC;
    PRB{TTI_counter}=prb_allocation;
  end

    rate_avg_collector=[rate_avg_collector;ue_rate_avg];
    [mue_out_sinr cue_out_sinr reue_out_sinr ue_out_sinr ue_out_efficiency mue_pos cue_pos reue_pos]=Figure_Throughput(UE_TP,UE_scheduled,UE_type,PRB,'user',antenna_configuration,1,UBS_layout);
    UE_efficiency_collector=[UE_efficiency_collector;ue_out_efficiency];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UE sinr%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UE_sinr_collector=[UE_sinr_collector;ue_out_sinr];
   
 
    
%%%%%%%%%%%%%%%%%%%%%%%UE pos - only scheduled UE%%%%%%%%%%%%%%%%%%%%%%%%%
    MUE_pos_collector=[MUE_pos_collector;mue_pos];
% end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Result%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
h1=cdfplot(UE_sinr_collector);
set(h1,'color','k');
hold on;


figure(3)
cdfplot(rate_avg_collector)
title('UE Throughput');


     



