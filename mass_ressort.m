

m = 1;
v0 = [-1 1];
k = 1;
dx = 0.01;
l = 2;

x = [-1 1];
r= 4;
h=plot(0,0,'MarkerSize',100,'Marker','.');
axis([-r r -r r]);
set(gca,'nextplot','replacechildren');



x(end+1,:) = x(end,:) - [1 -1]*k/(2*m)*dx^2*(x(end,1)-x(end,2)+l) + dx*v0 ;

set(h,'XData',x(end,:),'YData',[0,0]);
drawnow;


for t = 1:10000
    
    x(end+1,:) = 2*x(end,:) - [1 -1]*k/m*dx^2*(x(end,1)-x(end,2)+l) - x(end-1,:);
    
    set(h,'XData',x(end,:),'YData',[0,0]);
    drawnow;
    %pause(0.1)
end

plot(x)
