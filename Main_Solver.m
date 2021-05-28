clear all
clc
format long
type double

%-------------------------------------------------------------------------%
% Pressão em 'GPa', velocidades em 'm/s', fonte em 'm³/s', viscosidade    %
% em 'cP';                                                                %
%-------------------------------------------------------------------------%

global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge kmap ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime nflagface auxface gamma ...
    wells porousarea source O1 P1 T1 ve21 ve11 theta21 theta11 eps knownb ...
    neta1 esuel caso F q fract formethod timedimkey phflw rowposit ...
    regionfract region faces neigh state jointnodes malha weightDMP

%-------------------------------------------------------------------------%
% Escolha o exemplo:
    caso = 'Caso_3_Firozabadi_Cruzando';
% Escolha a malha:
    malha = 'Caso_3_Firozabadi_Cruzando'; 
%-------------------------------------------------------------------------%

[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,...
normals,esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
satlimit,pormap,bcflag,~,totaltime,kmap,wells,porousarea,ds,tconv,...
fract,elem_post,regionfract,coord_post,formethod,timedimkey,phflw,knownb,...
rowposit,region,numsline,substeps,maxsteps,benchmark,jointnodes] = preprocessor_frac2d;

%-------------------------------------------------------------------------%
 
%-------------------------------------------------------------------------%
eps = 1e-8; nw = 2; no = 2; state = {}; time = 0; timeref = 0; 
step = 0; VPI(1) = 0; cumulateoil(1) = 0; oilrecovery(1) = 1; 
finaltime = totaltime(2); dt_ref = 0; countime = 0; q = 0; watercut(1,1) = 0;
courant = 4*(ds==1)+0.9*(ds==2||ds==3); CFL = courant*(ds==1||ds==2)+4*(ds==3); 
bflux = zeros(size(bedge,1),1); influx = zeros(size(inedge,1),1);
if wells(1,1)==0 && ds==1, ds=2; end

[ source ] = sourceterm( elem );
nflag = calflag;
if ds==1, metsat='SEQ'; elseif ds==2, metsat='IMPES'; elseif ds==3, metsat='Streamline'; end 

% Outros Parâmetros
if strcmp(formethod,'NLFVPP')==1 || strcmp(formethod,'MPFAQL')==1 ...
   || strcmp(formethod,'NLFVDMP')==1 || strcmp(formethod,'MPFAH')==1
    [csi, p, tol, nit, nflagface, auxface, gamma] = preNLFV(kmap,nflag); 
    [weightDMP] = weightnlfvDMP(kmap);
else
    csi=0;tol=0;p=0;nit=0;weightDMP=0;nflagface=0;auxface=0;gamma=0;
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Cálculo dos parâmetros geométricos
[ O1, P1, T1, ve21, ve11, theta21, theta11, neta1 ] = geomparamLPEW2( coord, esurn2 );

% Parâmetros para cálculo de pressão
[Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn(kmap);

% Parâmetros para cálculo de saturação
[N,F,V,~,esuel,~,~,~,S_old,S_cont] = presaturation(wells,nflag);

% Preprocessor de Streamline
if ds==3
    [faces, neigh] = teste_unified_preprocess;
    if phflw==2, phflw = 1; end
end
%-------------------------------------------------------------------------%

% Verificação de Condição Inicial ----------------------------------------%
[ S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time ] = ...
  condinicial( S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time, metsat );
   
vpi_old = VPI(size(VPI,2)); v = 100*vpi_old; time2 = ceil(v)/100; pflag = 0;
%-------------------------------------------------------------------------%

%% IMPES / SEQ
while timeref < finaltime
%while t_old<totaltime  

    tic

    %% Cálculo da Pressão Implicta
    step = step + 1;
    
    % cálculo das mobilidades
    [mobility] = mobilityface(S_old,nw,no,S_cont);

    % cálculo da pressão
    [w, s, p, influx, bflux, q] = solverpressure(kmap, mobility,wells,weightDMP,...
                                  S_old, nflag, V, csi, N, nw, no, Hesq, Kde,...
                                  nit, tol, Kn, Kt, Ded, p, influx, bflux, q, step);

    %% Cálculo da Saturação
    if phflw~=0  
        
        % cálculo de passo de tempo (em 's'):------------------------------
        d_t = timestep(S_old,influx,bflux,CFL,nw,no,inedge,bedge);
        
        if step==1
            dt_ref = d_t;
            if strcmp(formethod,'NLFVDMP')==1||strcmp(formethod,'NLFVPP')==1
                d_t = 1e-12;
            end
        end  
        
        if step>1 && (d_t<(dt_ref/4)) && ds==1
            d_t = dt_ref/4; % Para evitar que o implícito dê passos muito pequenos.
        end
        %------------------------------------------------------------------
        
        % cálculo do fluxo fracional
        [f_elem] = fractionalflow(S_old,nw,no);

        % Saturação                    
        [ S_old ] = calcsaturation ( f_elem,S_cont,satlimit,visc,S_old,...
                    influx,bflux,d_t,wells,q,nw,no,elem,bedge,inedge,...
                    elemarea,pormap,ds,numsline,substeps,maxsteps,region );
                
        %% Reporte de Produção 
        [VPI,oilrecovery,cumulateoil,watercut,countime]=reportproduction(countime,...
            VPI,wells,f_elem,step,oilrecovery,cumulateoil,watercut,q,d_t,porousarea,S_old);

        vpi_old = VPI(step+1);
        countime_old = countime(step+1);
        
        if strcmp(timedimkey,'vpi')==1
            timeref = vpi_old;
        else
            timeref = countime_old;
        end
        
    else
       timeref = finaltime;
       d_t = 0;
    end
   
    fprintf ('- Step %u ------------------------------------\n',step);
    d_t
    
    satur = [ max(S_old) min(S_old) ]
    
    if max(S_old)>1.0000001 || min(S_old)<0
       supmax=find(S_old>1.0000001);
       S_old(supmax)
       elem(supmax,:)
       submin=find(S_old<0);
       S_old(submin)
       elem(submin,:)
       pause
    end 
    
    PVI = vpi_old
    
    %% Visualização
    if vpi_old >= time2
                   
        if ~exist(sprintf('%s\\%s\\%s\\%s\\Results',caso,malha,formethod,metsat),'dir')
             mkdir(sprintf('%s\\%s\\%s\\%s\\Results',caso,malha,formethod,metsat));
        end
            
        if max(S_old)<=1.0000001 && min(S_old)>=0
            
            save(sprintf('%s\\%s\\%s\\%s\\Results\\Saturation',caso,malha,formethod,metsat),'S_old');
            save(sprintf('%s\\%s\\%s\\%s\\Results\\VPI',caso,malha,formethod,metsat),'VPI');
            save(sprintf('%s\\%s\\%s\\%s\\Results\\Countime',caso,malha,formethod,metsat),'countime');
            save(sprintf('%s\\%s\\%s\\%s\\Results\\TimeStep',caso,malha,formethod,metsat),'step');
            save(sprintf('%s\\%s\\%s\\%s\\Results\\Watercut',caso,malha,formethod,metsat),'watercut');
            save(sprintf('%s\\%s\\%s\\%s\\Results\\OilRecovery',caso,malha,formethod,metsat),'oilrecovery');
            save(sprintf('%s\\%s\\%s\\%s\\Results\\CumulateOil',caso,malha,formethod,metsat),'cumulateoil');
            save(sprintf('%s\\%s\\%s\\%s\\Results\\TIME',caso,malha,formethod,metsat),'time');
            
            CopiaDeSeguranca( S_old,VPI,countime,step,watercut,oilrecovery,cumulateoil,time,metsat );
            
            if step==1
                save(sprintf('%s\\%s\\%s\\%s\\Results\\DT',caso,malha,formethod,metsat),'dt_ref');
            end
            
        end      
        
        passo = 1000*time2;
        postprocessor(p, S_old, passo, coord_post, elem_post,metsat)
        if ds==3
            vtpWriter(sprintf('%s\\%s\\%s\\%s\\Results\\PressOut',caso,malha,formethod,metsat),...
                passo, state.xyz, state.tof, state.time,'saturation', state.ssat);   
        end
        time2 = time2 + 0.01;

    end

    time = time+toc
    fprintf ('----------------------------------------------\n');
    
end
