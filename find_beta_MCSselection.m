

function [beta_LA]=find_beta_MCSselection(cqi_index,nue)

      beta_LA=zeros(nue,1);
      
  for ue=1:nue

      if cqi_index(ue)==1
           beta_LA(ue)=1.6;
      elseif cqi_index(ue)==2
          beta_LA(ue)=1.6;
      elseif cqi_index(ue)==3
          beta_LA(ue)=4.5;
      elseif cqi_index(ue)==4
          beta_LA(ue)=6.7;
      elseif cqi_index(ue)==5
          beta_LA(ue)=1.6;
      elseif cqi_index(ue)==6
          beta_LA(ue)=4.5;
      elseif cqi_index(ue)==7
          beta_LA(ue)=6.7;
      elseif cqi_index(ue)==8
          beta_LA(ue)=1.6;
      elseif cqi_index(ue)==9
          beta_LA(ue)=4.5;
      elseif cqi_index(ue)==10
          beta_LA(ue)=6.7;
      elseif cqi_index(ue)==11
          beta_LA(ue)=1.6;
      elseif cqi_index(ue)==12
          beta_LA(ue)=4.5;
      elseif cqi_index(ue)==13
          beta_LA(ue)=6.7;
      elseif cqi_index(ue)==14
          beta_LA(ue)=4.5;
      elseif cqi_index(ue)==15
          beta_LA(ue)=6.7;
      end
      
  end
  beta_LA;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
      