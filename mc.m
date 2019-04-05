% Calcul le centre de masse de la prot�ine
function [xs,ys] = mc(x,y,N,t)

xs = sum(x(t,:))/N;
ys = sum(y(t,:))/N;

end