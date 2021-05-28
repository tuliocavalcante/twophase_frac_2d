function [dfdS] = calcdfdS1(Sw,nw,no)
%Define global parameters
global satlimit visc

%Define "dfdS" and "dgamadS" according to "dertype"
%Initialize some propoerties (two-phase flow)
Swi = satlimit(1);
Sor = satlimit(2);
miw = visc(1);
mio = visc(2);
%Initialize "dfdS" and "dgamadS"
dfdS = zeros(length(Sw),1);

%Fit parameter (water and oil)
kwmax = 1;
komax = 1;

%Calculate the derivate
for i = 1:length(Sw)
    %Define some terms:
    term1 = 1 - Swi - Sor;
    term2 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1));
    term3 = komax*((1 - ((Sw(i) - Swi)/term1))^no)/mio;
    term4 = kwmax*(((Sw(i) - Swi)/term1)^nw)/miw;
    term5 = kwmax*(((Sw(i) - Swi)/term1)^nw);
    term6 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1))/(term1*miw);
    term7 = no*komax*((1 - ((Sw(i) - Swi)/term1))^(no - 1))/...
        (term1*mio);
    %Calculate "dfdS"
    dfdS(i) = (term2/(term1*(term3 + term4)*miw)) - ...
        (term5*(term6 - term7))/(((term3 + term4)^2)*miw);
end  %End of FOR

