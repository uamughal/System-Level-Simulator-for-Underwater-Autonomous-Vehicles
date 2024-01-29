% function [mbs_layout area_grid grid_cell_ind] = MBS_Layout( isd,n_tier,tower_height,isplot,sim_res,mue_free )

function [UBSC_layout area_grid grid_cell_ind]=UBSC_Layout(ISD,n_tier,isplot,sim_res,UBS_free);
r=ISD/2;
UBSC_coordinate(1:1)=0;
n_cells=1;

UBSC_layout=zeros(n_cells,8);
UBSC_layout(:,1)=real(UBSC_coordinate);
UBSC_layout(:,2)=imag(UBSC_coordinate);


% generate area grid
ag_x=-r:sim_res:r;
ag_y=-r:sim_res:r;
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
area_grid=repmat(ag_pos,n_cells,1)+reshape(repmat(UBSC_coordinate(1:3:end),length(ag_pos),1),n_cells*length(ag_pos),1)+r/2;
grid_cell_ind=[cell_ind;cell_ind+3;cell_ind+6;cell_ind+9;cell_ind+12;cell_ind+15;cell_ind+18];

% plot hexagonal grid for illustration
if isplot==1
    X=[];
    Y=[];
    Ax=[];
    Ay=[];
    figure
    hold on
    grid on
  
%   X_base=[-r,-r/2,r/2,r,r/2,-r/2;
%         -r/2,r/2,r,r/2,-r/2,-r];
%     Y_base=[0,   r, r, 0,    -r, -r;
%        r, r, 0,   -r, -r, 0];

%    
%     Ax_base=[0,     0,       0;
%              1,   -1/2,     -1/2]*r;
%     Ay_base=[0,     0,        0;
%              0,     1,-1]*r;
    Ax_base=[0,     0,       0;
             1,   -1/2,     -1/2]*r/3;
    Ay_base=[0,     0,        0;
             0,     sqrt(3)/2,-sqrt(3)/2]*r/3;
 
%     Ax_base=[0,     0,       0;
%              0,   0,     0];
%     Ay_base=[0,     0,        0;
%              0,     0,0];
%          
   for ii=1:n_cells
%         X=[X,X_base];
%         Y=[Y,Y_base];
%         X=[X,X_base+real(BS_coordinate(ii))-1/2*r/2];
%         Y=[Y,Y_base+imag(BS_coordinate(ii))+sqrt(3)/2*r/2];
%            X=[X,X_base+real(BS_coordinate(ii))-1/2*r/2];
%            Y=[Y,Y_base+1];
%         X=[X,X_base+real(BS_coordinate(ii))-1/2*r/2];
%         Y=[Y,Y_base+imag(BS_coordinate(ii))-sqrt(3)/2*r/2];
        Ax=[Ax,Ax_base+real(UBSC_coordinate(ii))];
        Ay=[Ay,Ay_base+imag(UBSC_coordinate(ii))];
   end
N_sides = 6;
t=(1/(N_sides*2):1/N_sides:1)'*2*pi;
x=sin(t);
y=cos(t);
X=r*[x; x(1)];
Y=r*[y; y(1)];
plot(x,y)
axis square
   line(X,Y,'Color','k','LineWidth',2);
   line(Ax,Ay,'Color','b','LineWidth',2);
   plot(real(UBSC_coordinate),imag(UBSC_coordinate),'bs','MarkerFaceColor','b');

end
