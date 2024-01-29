function pl=Path_Loss(UBSC_layout,UBS_layout)
% intialize distance matrix
dist=zeros(size(UBS_layout,1),size(UBSC_layout,1));     
% intialize PL matrix
pl=dist; 
% pl_a=dist;
pl_t=dist;        

% number of BSs and UEs
nbs=size(UBSC_layout,1);
nue=size(UBS_layout,1);
dist_x=repmat(UBS_layout(:,1),1,nbs)-repmat(UBSC_layout(:,1).',nue,1);
dist_y=repmat(UBS_layout(:,2),1,nbs)-repmat(UBSC_layout(:,2).',nue,1);
dist_z=repmat(UBS_layout(:,3),1,nbs)-repmat(UBSC_layout(:,3).',nue,1);
dist=sqrt(dist_x.^2+dist_y.^2+dist_z.^2);
dist_a=log10(dist);

% find the UE-BS matrix index for each individual propagation type
ue_ind=find(UBS_layout(:,3)==0);
bs_ind=find(UBSC_layout(:,3)==0);

n_ue=length(ue_ind);
n_bs=length(bs_ind);
if n_ue~=0 && n_bs~=0
    ue_ind=repmat(ue_ind,n_bs,1);
    bs_ind=reshape(repmat(bs_ind.',n_ue,1),n_ue*n_bs,1);
    index=sub2ind([nue nbs],ue_ind,bs_ind);
    pl_t(index)=0; 
end


%for pico&macro scenario
ue_ind=find(UBS_layout(:,6)==0);

pibs_ind=find(UBSC_layout(:,4)~=0);
mbs_ind=find(UBSC_layout(:,4)==0);
n_pibs=length(pibs_ind);
n_mbs=length(mbs_ind);
n_ue=length(ue_ind);
if n_mbs~=0
    
   
    
    % For UBSC and UBS
    % UBS
    ue_ind=repmat(ue_ind,n_pibs,1);
    pibs_ind=reshape(repmat(pibs_ind.',n_ue,1),n_ue*n_pibs,1);
    index=sub2ind([nue nbs],ue_ind,pibs_ind);
    pl_t(index)=5;
    % UBSC
    ue_ind=find(UBS_layout(:,6)==0);
    n_ue=length(ue_ind);
    ue_ind=repmat(ue_ind,n_mbs,1);
    mbs_ind=reshape(repmat(mbs_ind.',n_ue,1),n_ue*n_mbs,1);
    index=sub2ind([nue nbs],ue_ind,mbs_ind);
    pl_t(index)=6;
    
    
end


%% for UBSC and UBS 
ind =find(pl_t==6);
if ~isempty(ind)
fkz=5.5*10^3;
xa = fkz.^2;
fa = (((0.11*(xa))./(1+(xa))) + ((44*(xa))./(4100+(xa))) + (2.75*(xa)/10000) + 0.003);
k=1.5;     % Spreading Factor
%  A_C=log10(fa);
  pl(ind)=((k*log10(dist_a(ind)))+(dist_a(ind)*log10(fa)));
%   pl_a(ind)=log10(pl(ind));
end
% 
% 
% ind =find(pl_t==5);
% if ~isempty(ind)
%     f=5.5*10^3;
%     Ab_Co=((0.11*f^2/(1+f^2))+(44*f^2/(1400+f^2))+((2.75*10^-4)*f^2)+0.003);   % Absorption Coeficient
%     A_C=log10(Ab_Co);
%     k=1.5;     % Spreading Factor
%     pl(ind)=(k*log10(dist(ind)))+(A_C*(dist(ind)));
%     pl_a(ind)=log10(pl(ind));
% end
% ind=find(pl_t==6);
% if ~isempty(ind)
%        f=5.5*10^3;
%     Ab_Co=((0.11*f^2/(1+f^2))+(44*f^2/(1400+f^2))+((2.75*10^-4)*f^2)+0.003);   % Absorption Coeficient
%     A_C=log10(Ab_Co);
%     k=1.5;     % Spreading Factor
%     pl(ind)=(k*log10(dist(ind)))+(A_C*(dist(ind)));
%     pl_a(ind)=log10(pl(ind));
% %     pl(ind)=(k*log10(dist(ind)))+(log10*Ab_Co(dist(ind)));
% end

    
