clc;clear all; close all;
%%

xi = 0.9;

k = 1;
h = 3;
dL = 5/100;
L = 30/100;
cv = 1000;
rho = 2000;
q = 2000;
Ta = -10;
Ti = 20;

dx = 0.01;
dt = 100*dx^2*cv*rho/k;

%% initialisation

c1=-k; 
c2=h; 
c3=-h*Ta;

d1=-k; 
d2=-h; 
d3=h*Ti;

N = L/dx;
xl=linspace(0, L , N)';

Sx = q./(1+((xl-0.5*L)/dL).^2);   

A = zeros(N);

v = [];

x = Ta + (Ti-Ta)/((2*k)/(h*L)+1)*(k/(h*L)+xl/L);


Tmax = zeros(1,N);
Tmax_0 = max(x);

%% Matrice A et b

j = 1;

for i = 2:N-1
    A(i,j) = 1 ;
    A(i,j+1) = -2 ;
    A(i,j+2) = 1;
    j = j + 1;
end
    b = -(Sx/k)  *dx^2;
    
%% Conditions frontières

A(1,1)      =  (2*c2*dx-3*c1);
A(1,2)      =  4*c1;
A(1,3)      =   -c1;
b(1)        =  -2*c3*dx;

A(N,N-2)    =  d1;
A(N,N-1)    =  -4*d1;
A(N,N)      =  3*d1+2*d2*dx;
b(N)        =  -2*d3*dx;


%% 
di = ones(1,N);
di(1) = 0; 
di(end) = 0;
M = diag(di);


%%
j = 1;

for i = 1:dt:1e6

bp  =  (M+( (dt*k)/(cv*rho*dx^2) ) * (1-xi)*A )*x-(dt*k)/(cv*rho*dx^2)*b;

Ap =(M - ( (dt*k)/(cv*rho*dx^2) ) * xi *A );

%% Résolution

x =  Ap\bp;

Tmax(j) = Tmax_0 +0.99*(max(x)-Tmax_0);
j = j+1;

end
%% Graph
t = 1:dt:1e6;
y = ones(1,length(t))*59.8622;
hold on
plot(t,Tmax)
plot(t,y)
hold off