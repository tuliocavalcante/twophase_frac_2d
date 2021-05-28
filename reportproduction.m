function [VPI,oilrecovery,cumulateoil,watercut,countime]=reportproduction...
         (countime,VPI,wells,f_elem,cont,oilrecovery,cumulateoil,watercut,...
         q,dt,porousarea,S_old)

if isempty(wells)==0 && max(wells(:,1))~=0
    
    sumvpi = 0;
    sumcumoil = 0;
    sumvol  = 0;
    sumwater=0;

    for iwell = 1:size(wells,1)
        
        if wells(iwell,3) == 0 % no poço produtor
            
            sumvol = sumvol + porousarea(wells(iwell,1));

            sumwater = sumwater + S_old(wells(iwell,1))*porousarea(wells(iwell,1));
                        
            f_o = 1 - f_elem(wells(iwell,1));
            
            sumcumoil = sumcumoil - f_o*q(wells(iwell,1))*dt;

        elseif wells(iwell,3) ~= 0 % no poço injetor
            
            sumvpi = sumvpi + f_elem(wells(iwell,1))*q(wells(iwell,1))*dt*(1/sum(porousarea));

        end
        
    end
    
    countime(cont+1) = countime(cont) + dt;
    
    VPI(cont+1) = VPI(cont) + sumvpi;

    cumulateoil(cont+1) = cumulateoil(cont) + sumcumoil;

    oilrecovery(cont+1) = 1 - (sumwater/sumvol);

    watercut(cont+1)= sumwater/sumvol;
            
elseif isempty(wells)~=0
    
    countime(cont+1) = countime(cont) + dt;
    
    VPI(cont+1) = VPI(cont) + sum((S_old.*porousarea')*(1/sum(porousarea)))*(1/size(S_old,1));

    cumulateoil(cont+1) = cumulateoil(cont) ;

    oilrecovery(cont+1) = 1 ;

    watercut(cont+1)= 0;

end