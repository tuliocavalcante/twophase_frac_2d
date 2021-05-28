function [f_elem] = fractionalflow(S_old,nw,no)

global visc satlimit phflw

if phflw==2
 
f_elem = zeros(size(S_old,1),1);

% loop de fases internas para calcular as mobilidade e fluxo fracional
for i = 1:size(S_old,1) 
    % calculo das permeabilidade usando modelo Brooks e Corey
     
    % calculando as permeabilidades relativas e a mobilidade do elemento

    Krw2 = ((S_old(i) - satlimit(1))/(1-satlimit(1)-satlimit(2)))^nw;
    
    Kro2 = ((1 - S_old(i) - satlimit(1))/(1-satlimit(1)-satlimit(2)))^no;
    
    L2 = Krw2/visc(1) + Kro2/visc(2);
         
    f_elem(i) = (Krw2/visc(1))/L2; 
 
end

elseif phflw==1
    
    f_elem = S_old;
    
end
       
end
