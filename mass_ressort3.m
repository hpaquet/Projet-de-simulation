clear all;close all; clc;
%%
m = 1;
v0 = [-1 1];
k = 1;
dt = 0.1;
l = 2;
g = 0;
xi = 1;
D = 1;


x = [-1 1];
r= 3;
h=plot(0,0,'MarkerSize',100,'Marker','.');
axis([-r r -r r]);
set(gca,'nextplot','replacechildren');

%%
% x(end+1,:) =  2*dx*force(x(end,1),x(end,2),k,l) - dx*v0 ;

C = 1/(4*m);
C1 = C*4*m;
C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(end,1),x(end,2),k,l));
C2 = C*(2*m+g*dt)*2*dt;
x(end+1,:) = C1.*x(end,:) + C3 + C2.*v0 ;

set(h,'XData',x(end,:),'YData',[0,0]);
drawnow;
%%

for t = 1:10000
    
%     x(end+1,:) = 2*force(x(end,1),x(end,2),k,l)*dx - x(end-1,:);
    C = 1/(2*m-g*dt);
    C1 = C*4*m;
    C2 = -C*(2*m+g*dt);
    C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(end,1),x(end,2),k,l));

    x(end+1,:) = C1*x(end,:) + C2*x(end-1,:) + C3;
    
    set(h,'XData',x(end,:),'YData',[0,0]);
    pause(0.1)
    drawnow;

end



%%

plot(x)

%%

function F = force(x1,x2,k,a)
    F = - [1 -1]*k*(x1-x2+a);
end