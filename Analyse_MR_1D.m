clc;close all;clear all;
%%
% Fonctionne parfaitement
%
%
%
%%

m = 1;
k = 10;
d = 0.1;
l = 4;
g = 0;
xi = 0;
D = 0;
tf = 100;
j = 1;

x_sol = [];

%%

for i = [1 0.5 0.25 0.125 0.0625 0.0312]
    
    dt = i*d;
    N = tf/i;
    
    for n = 1:2
        x(1,n) = n;
        v0(1,n) = (-1)^(n-1);
    end
    
    for n = 1:2
        C = 1/(4*m);
        C1 = C*4*m;
        C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(1,:),n,k,l));
        C2 = C*(2*m+g*dt)*2*dt;
        x(2,n) = C1.*x(1,n) + C3 + C2.*v0(n) ;
    end
    
    
    for t = 2:N
        
        for n =1:2
            C = 1/(2*m-g*dt);
            C1 = C*4*m;
            C2 = -C*(2*m+g*dt);
            C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(t,:),n,k,l));
            
            x(t+1,n) = C1*x(t,n) + C2*x(t-1,n) + C3;
            
        end
        
    end
    
    x_sol(end+1,:) = max( x - x(1,:));
    
end

err = [];

for k = 1:length(x_sol)-1
    err(k,:) = abs(x_sol(k,:) - x_sol(k+1,:));
end

N = [0.2000,0.4000,0.8000,1.6000,3.2051]'*1e3;
% 
loglog(N,err)


%%

function F = force(x,n,k,a)

F = 0;

if n-1 >= 1
    F = F - k*(x(n)-x(n-1)-a);
end
if n+1 <= length(x)
    F = F - k*(x(n)-x(n+1)+a);
end

end
