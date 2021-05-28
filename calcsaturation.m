function [ S_old ] = calcsaturation ( f_elem,S_cont,satlimit,visc,S_old,...
                     influx,bflux,d_t,wells,q,nw,no,elem,bedge,inedge,...
                     elemarea,pormap,ds,numsline,substeps,maxsteps,region )
%
       
if ds==2
    
    [S_old]= firstorderstandard(S_old,influx,bflux,q,f_elem,d_t,wells,...
                                S_cont,nw,no,elem,inedge,bedge,pormap,...
                                elemarea,satlimit,visc,region);             % IMPES
      
elseif ds==1
    
    [S_old] = firstorderstdseqimp(S_old,influx,bflux,d_t,wells,q,nw,no,...
                                 elem,bedge,inedge,elemarea,pormap,region); % SEQUENCIAL IMPLÍCITO
      
elseif ds==3
    
    [S_old] = firstorderstreamline(influx,bflux,S_old,numsline,substeps,maxsteps,d_t,nw,no); % STEAMLINES
    
end

            
end