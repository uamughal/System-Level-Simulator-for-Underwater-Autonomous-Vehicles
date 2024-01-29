
function pl_slow=fd_ch(UBSC_layout,UBS_layout)
ue_ind=find(UBS_layout(:,3)==0);
bs_ind=find(UBSC_layout(:,3)==0);

ts = 1/5000;    % sampling period
fd =4;
UBS_layout(:,8)=[];
UBSC_layout(:,8)=[];
values=UBS_layout(ue_ind,:);UBSC_layout(bs_ind,:);
dist11= reshape(values,1,147);      
[ch]=gen_ch(ts,fd);    
pl_slow = filter( ch,dist11);
end