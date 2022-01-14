function [c ceq] = chooseGoatDaysConstrain(G)

for i=1:20
    c(i) = -G(i);
    c(i+20) = G(i)-1;
end
ceq=[];