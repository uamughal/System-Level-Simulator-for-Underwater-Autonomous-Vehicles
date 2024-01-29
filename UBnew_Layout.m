function [UBSC_layout area_grid grid_cell_ind]=UBnew_Layout(isd,n_tier,isplot,sim_res,UBS_free);
%  r=isd/2;
% r=2/3*isd;
r=0.575*isd;
UBSC_coordinate(1)=0;
n_cells=1;

if n_tier>=1
% first ring
n_cells=7;
UBSC_coordinate(2)=UBSC_coordinate(1)+(isd*sqrt(3)/2)+(isd*1/2*1j);
UBSC_coordinate(3)=UBSC_coordinate(1)+isd*1j;
UBSC_coordinate(4)=UBSC_coordinate(1)-isd*sqrt(3)/2+1j*isd*1/2;
UBSC_coordinate(5)=UBSC_coordinate(1)-isd*sqrt(3)/2-1j*isd*1/2;
UBSC_coordinate(6)=UBSC_coordinate(1)-1j*isd;
UBSC_coordinate(7)=UBSC_coordinate(1)+isd*sqrt(3)/2-1j*isd*1/2;
end



UBSC_layout=zeros(n_cells,8);
UBSC_layout(:,1)=real(UBSC_coordinate);
UBSC_layout(:,2)=imag(UBSC_coordinate);


% generate area grid
% ag_x=-r:sim_res:r;
% ag_y=-r:sim_res:r;
ag_x=-r/2:sim_res:r/2;
ag_y=-r*sqrt(3)/4:sim_res:r*sqrt(3)/4;
N_x=length(ag_x);
N_y=length(ag_y);
ag_pos=[reshape(repmat(ag_x,N_y,1),N_x*N_y,1),repmat(ag_y',N_x,1)];
ind=(ag_pos(:,2)-sqrt(3)*(ag_pos(:,1)-r/2)>0)&...
    (ag_pos(:,2)-sqrt(3)*(ag_pos(:,1)+r/2)<0)&(ag_pos(:,2)-sqrt(3)*...
    (r/2-ag_pos(:,1))<0)&(ag_pos(:,2)+sqrt(3)*(ag_pos(:,1)+r/2)>0);
ag_pos=ag_pos(ind,:);
ag_pos=[ag_pos;ag_pos(:,1)-3/4*r,ag_pos(:,2)+sqrt(3)/4*r;...
    ag_pos(:,1)-3/4*r,ag_pos(:,2)-sqrt(3)/4*r];
cell_ind=reshape(repmat([1 2 3],size(ag_pos,1)/3,1),size(ag_pos,1),1);
ind= abs(ag_pos(:,1)+1j*ag_pos(:,2)+r/2)>=UBS_free;
ag_pos=ag_pos(ind,:);
cell_ind=cell_ind(ind);
ag_pos=ag_pos(:,1)+1j*ag_pos(:,2);
area_grid=repmat(ag_pos,n_cells,1)+reshape(repmat(UBSC_coordinate(1:1:end),length(ag_pos),1),n_cells*length(ag_pos),1)+r/2;
grid_cell_ind=[cell_ind;cell_ind+1;cell_ind+2;cell_ind+3;cell_ind+4;cell_ind+5;cell_ind+6];

% plot hexagonal grid for illustration
if isplot==1
    X=[];
    Y=[];
    Ax=[];
    Ay=[];
    figure
    hold on
    grid on
%   
%   X_base=[-r,-r/2,r/2,r,r/2,-r/2;
%         -r/2,r/2,r,r/2,-r/2,-r];
%     Y_base=[0,   r, r, 0,    -r, -r;
%        r, r, 0,   -r, -r, 0]*sqrt(3)/2;


    Ax_base=[0,     0,       0;
             1,   -1/2,     -1/2]*r;
    Ay_base=[0,     0,        0;
             0,     sqrt(3)/2,-sqrt(3)/2]*r;
         
 N_sides = 6;
t=(1/(N_sides*2):1/N_sides:1)'*2*pi;  %length of hexagonal lines
x=sin(t);
y=cos(t);
X_base=(r)*[x; x];      %we have to only consider r for further
                        %changes in coding. Leave the rest of the 
                        %things as it is.
Y_base=(r)*[y; y];
          
   for ii=1:n_cells
%         X=[X,X_base];
%         Y=[Y,Y_base];
%         X=[X,X_base+real(BS_coordinate(ii))-1/2*r/2];
%         Y=[Y,Y_base+imag(BS_coordinate(ii))+sqrt(3)/2*r/2];
%            X=[X,X_base+real(BS_coordinate(ii))-1/2*r/2];
%            Y=[Y,Y_base+1];
%         X=[X,X_base+real(BS_coordinate(ii))-1/2*r/2];
%         Y=[Y,Y_base+imag(BS_coordinate(ii))-sqrt(3)/2*r/2];
        X=[X,X_base+real(UBSC_coordinate(ii))];
        Y=[Y,Y_base+imag(UBSC_coordinate(ii))];
        Ax=[Ax,Ax_base+real(UBSC_coordinate(ii))];
        Ay=[Ay,Ay_base+imag(UBSC_coordinate(ii))];
   end
   
% axis square
   line(X,Y,'Color','k','LineWidth',3);
   line(Ax,Ay,'Color','c','LineWidth',0.1);
   plot(real(UBSC_coordinate),imag(UBSC_coordinate),'bs','MarkerFaceColor','b');

end
