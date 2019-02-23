clear all;close all; clc;
%%
m = 1;

k = 5;
dt = 0.1;
l = 2;
g = 0;
xi = 1;
D = 1;

N=64;

for n = 1:N
    x(1,n) = 2*(n);
    v0(1,n) = (-1)^n;
end
% x = [-1 1];
r= 3;
h=plot(0,0,'MarkerSize',100,'Marker','.');
axis([-0 100 -1 1]);
set(gca,'nextplot','replacechildren');

%%
% x(end+1,:) =  2*dx*force(x(end,1),x(end,2),k,l) - dx*v0 ;
for n = 1:N
    C = 1/(4*m);
    C1 = C*4*m;
    C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(1,:),n,k,l));
    C2 = C*(2*m+g*dt)*2*dt;
    x(2,n) = C1.*x(1,n) + C3 + C2.*v0(n) ;
end

set(h,'XData',x(end,:),'YData',zeros(1,N));
drawnow;
%%

for t = 2:10000
    
    %     x(end+1,:) = 2*force(x(end,1),x(end,2),k,l)*dx - x(end-1,:);
    for n =1:N
        C = 1/(2*m-g*dt);
        C1 = C*4*m;
        C2 = -C*(2*m+g*dt);
        C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(t,:),n,k,l));
        
        x(t+1,n) = C1*x(t,n) + C2*x(t-1,n) + C3;
        
    end
    
    set(h,'XData',x(end,:),'YData',zeros(1,N));
    pause(0.1)
    drawnow;
    
end



%%

plot(x)

%%

function F = force(x,n,k,a)
    
    F = 0;
    
    if n-1 >= 1
        F = F + k*(x(n-1)-x(n)+a);
    end
    if n+1 <= length(x)
        F = F - k*(x(n)-x(n+1)+a);
    end
    
end
