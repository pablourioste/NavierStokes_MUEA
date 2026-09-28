function [Ru,Rv]=computeR(V,cu,cv,du,dv)

Ru=-cu/V+du/V;
Rv=-cv/V+dv/V;

end