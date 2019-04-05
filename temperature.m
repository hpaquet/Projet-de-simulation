% Calcule la température du système
function T = temperature(x,y,N,dt,m,t)

C = m/(2*physconst('Boltzmann')*N);

vx = 3*x(t,:)-4*x(t-1,:)+x(t-2,:);
vy = 3*y(t,:)-4*y(t-1,:)+y(t-2,:);
v2 = vx.^2+vy.^2;

T = C*sum(v2);

end