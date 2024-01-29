function [ue_layout,ap_layout,UE_AP]=BS_Association(UE_Layout,MBS_layout)

% UBS association
N_UE=size(UE_Layout,1);
N_BS=size(MBS_layout,1);
UE_AP=zeros(N_UE,N_BS,'int8');
% switch_ue=[];
% the indoor-ue associate with MBS when there is no HBS
% if strcmp(Indoor_Switch_Type,'MBS')
%     for ii=1:size(UE_Layout,1)
% %         block_ind=UE_Layout(ii,8);% building index 1*1
% %         bs_ind=find(AP_Layout(:,8)==block_ind);%bs_indexes in the building
%         % check the flat number of ue, compare with those in HBS
%         bs_served=find(UE_Layout(ii,7)==AP_Layout(bs_ind,7), 1);%the hbs in same apartment with ue
%         if isempty(bs_served)
%             switch_ue=[switch_ue ii]; % indexes of indoor UEs with no hbs in the same apartment
%         else
%             UE_AP(ii,bs_ind(bs_served))=1;
%         end
%     end
% end
% UE_Layout(switch_ue,7)=0; %set the apartments indexes 0, for there is no hbses inside the apartment to serve those indoor UEs

% add marocell UE
if isempty(UE_Layout)
    ue_layout=MC_UE;
else
    ue_layout=[UE_Layout];
end

% add macrocell BS
ap_layout=MBS_layout;