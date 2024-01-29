% function [UBS_layout]=UBS_Layout(UBSC_layout,isplot,NUE,dropf,area_grid,cell_ind,N_UBS,n_cells,ISD,D1,D2,D3)
function [UBS_layout]=UBS_Layout(UBSC_layout,D1,D2,D3)

cell_pos=UBSC_layout(:,1)+1j*UBSC_layout(:,2);
nUBS_cell=length(cell_pos);
UBS_layout=zeros(nUBS_cell*3,8);
d_UBS_pos=[D2, -D1*(1/2)+1j*D1*sqrt(3)/2, -D3*(1/2)-1j*D3*sqrt(3)/2];
UBS_pos_a=repmat(d_UBS_pos,nUBS_cell,1)+repmat(cell_pos,1,3);
UBS_pos=reshape(UBS_pos_a,21,1);
% UBS_pos=reshape(UBS_pos_a,12,1);
UBS_layout(:,1)=real(UBS_pos);
UBS_layout(:,2)=imag(UBS_pos);
figure(1)
plot(UBS_layout(:,1),UBS_layout(:,2),'r*')
hold on
