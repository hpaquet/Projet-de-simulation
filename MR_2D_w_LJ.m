clc;close all;clear all;
%%

global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);

%%

m = 1;

k = 1;

kb = 1.38e-23;
dt = 0.001;
l = 1;
g = 1e-6;
xi = 1;
D = 1/g;

N=3;

i = 0;

T = [];

for n = 1:N
    x(1,n) = randi(10);
    y(1,n) = randi(10);
    v0(1,n) = (-1)^(randi(N));
    v0(2,n) = (-1)^(randi(N));
end

xx = 20;
yy = 20;

figure(1);
f1 = subplot(1,1,1);
h1=plot(0,0,'MarkerSize',100,'Marker','.','LineWidth',5);
axis(f1,[-xx xx -yy yy]);
set(gca,'nextplot','replacechildren');
figure(2)
h2=animatedline;


%%

for n = 1:N
    
    xi = normrnd(0,1);
    
    C = 1/(4*m);
    C1 = C*4*m;
    C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(1,:),y(1,:),n,l));
    C2 = C*(2*m+g*dt)*2*dt;
    x(2,n) = round(C1.*x(1,n) + C3(1) + C2.*v0(1,n),3);
    y(2,n) = round(C1.*y(1,n) + C3(2) + C2.*v0(2,n),3);
    
end

set(h1,'XData',x(end,:),'YData',y(end,:));

drawnow;

%%

for t = 2:10000
    
    for n =1:N
        
        xi = normrnd(0,1);
        
        C = 1/(2*m-g*dt);
        C1 = C*4*m;
        C2 = -C*(2*m+g*dt);
        C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(t,:),y(t,:),n,l));
        
        x(t+1,n) = round(C1.*x(t,n) + C3(1) + C2*x(t-1,n),3);
        y(t+1,n) = round(C1.*y(t,n) + C3(2) + C2*y(t-1,n),3);
        
    end
    
    T = temperature(x,y,N,kb,dt,m);
    [xs,ys] = mc(x,y,N);
    
    if mod(i,50) == 0
        addpoints(h2,t,T);

        set(h1,'XData',x(end,:),'YData',y(end,:));
        axis(f1,[-xx+xs xx+xs -yy+ys yy+ys]);

        drawnow;
    
    end
    
    if KEY_IS_PRESSED
        close all;
        break;
    end
    
end


%%

function [xs,ys] = mc(x,y,N)
    xs = sum(x(end,:))/N;
    ys = sum(y(end,:))/N;
end


function F = force(x,y,n,a)

F = [0 0];

for i = 1:length(x)
    if i == n-1
        dx = x(i)-x(n);
        dy = y(i)-y(n);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        k = LenardJones(a,delta,10,2);
        
        F(1) = F(1) + k*((delta)*cos(theta));
        F(2) = F(2) + k*((delta)*sin(theta));
        
    elseif i == n+1
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        k = LenardJones(a,delta,10,2);
        
        F(1) = F(1) - k*((delta)*cos(theta));
        F(2) = F(2) - k*((delta)*sin(theta));
        
    elseif i ~= n
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        f = LenardJones(a,delta,-0.3,1);
        
        F(1) = F(1) - f*cos(theta);
        F(2) = F(2) - f*sin(theta);
        
    end
end
end


function T = temperature(x,y,N,k,dt,m)

T = 0;
C = m/(8*k*N*dt^2);

for i = 1:N
    
    vx = 3*x(end,i)-4*x(end-1,i)+x(end-2,i);
    vy = 3*y(end,i)-4*y(end-1,i)+y(end-2,i);
    
    v2 = vx^2+vy^2;
    
    T = T + v2;
    
end

T = C*T;

end


function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
end
