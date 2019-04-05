clc;close all;clear all;
%%


global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);

%%

m = 1;
kb = 1.38e-23;
k = 1;
T = 300;
dt = 0.01;
l = 3;

N= 3;

for n = 1:N
    x(1,n) = 2*n;
    y(1,n) = 0;%2*(n-1);
    v0(1,n) = 0;%(-1)^(n);
    v0(2,n) = 0;
end

subplot(2,2,1:2)
h=plot(0,0,'MarkerSize',100,'Marker','.','LineWidth',5);
axis([-N 4*N -5 10]);
set(gca,'nextplot','replacechildren');
subplot(2,2,3)
for i =1:N
    pl1(i) = animatedline('color',rand(1,3));
end
axis([0 100 -0 10]);
subplot(2,2,4)
for i =1:N
    pl2(i) = animatedline('color',rand(1,3));
end
axis([0 100 -5 10]);
%%

for n = 1:N
    C3 = dt^2*force(x(1,:),y(1,:),n,l)/(2*m);

    x(2,n) = x(1,n) + C3(1) + dt*v0(1,n);
    y(2,n) = y(1,n) + C3(2) + dt*v0(2,n);
end

set(h,'XData',x(end,:),'YData',y(end,:));

drawnow;

%%

for t = 2:10000
    
    for n =1:N
        
        C3 = dt^2*force(x(t,:),y(t,:),n,l)/m;
        
        x(t+1,n) = 2*x(t,n) + C3(1) - x(t-1,n);
        y(t+1,n) = 2*y(t,n) + C3(2) - y(t-1,n);
        
    end
    
    set(h,'XData',x(end,:),'YData',y(end,:));
    
    for i = 1:N
        addpoints(pl1(i),t,y(end,i));
        axis([t-20 20+t -0 10]);
    end
    
    
    for i = 1:N
        addpoints(pl2(i),t,x(end,i));
        axis([t-20 20+t -5 10]);
    end
    
    
    %pause(0.1)
    drawnow;
    
    if KEY_IS_PRESSED
        close all;
        break;
    end
    
end


%%

function F = force(x,y,n,a)

F = [0 0];
for i = 1:length(x)
    if i == n-1
        dx = x(i)-x(n);
        dy = y(i)-y(n);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        k = LenardJones(a,delta,1,2)*10;
        
        F(1) = F(1) + k*((delta)*cos(theta));
        F(2) = F(2) + k*((delta)*sin(theta));
        
    elseif i == n+1
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        k = LenardJones(a,delta,1,2)*10;
        
        F(1) = F(1) - k*((delta)*cos(theta));
        F(2) = F(2) - k*((delta)*sin(theta));
        
    elseif i ~= n
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        f = LenardJones(a,delta,1,1);
        
        F(1) = 0;%F(1) - f*cos(theta);
        F(2) = 0;%F(2) - f*sin(theta);
    end
end
end
%%

function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
end
