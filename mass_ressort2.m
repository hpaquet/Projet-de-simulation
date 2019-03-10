clear all;close all; clc;
%%
m = 1;
v0 = [-1 1];
k = 1;
dt = 0.1;
l = 2;

x = [-1 1];
r= 4;
h=plot(0,0,'MarkerSize',100,'Marker','.');
axis([-r r -r r]);
set(gca,'nextplot','replacechildren');

%%
% x(end+1,:) =  2*dx*force(x(end,1),x(end,2),k,l) - dx*v0 ;

C1 = 1;
C3 = force(x(end,1),x(end,2),k,l)*dt^2/m;
C2 = dt;
x(end+1,:) = C1*x(end,:) + C3 + C2*v0 ;

set(h,'XData',x(end,:),'YData',[0,0]);
drawnow;
%%

for t = 1:10000
    
%     x(end+1,:) = 2*force(x(end,1),x(end,2),k,l)*dx - x(end-1,:);

    C1 = 2;
    C2 = -1;
    C3 = force(x(end,1),x(end,2),k,l)*dt^2/m;

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
