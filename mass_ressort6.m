clc;close all;clear all;
%%

global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);

%%

m = 10;

k = 10;

dt = 0.1;
l = 2;
g = 0;
xi = 0;
D = 0;

N=64;

for n = 1:N
    x(1,n) = 2*(n);
    y(1,n) = 0;
    v0(1,n) = (-1)^n;
    v0(2,n) = (-1)^(n-1);
end

% subplot(2,1,1)
h=plot(0,0,'MarkerSize',100,'Marker','.','LineWidth',5);
axis([-10 140 -10 10]);
set(gca,'nextplot','replacechildren');
% subplot(2,1,2)
% for i =1:N
%     pl(i) = animatedline('color',rand(1,3));
% end
% axis([0 100 -0 10]);
%%

for n = 1:N
    C = 1/(4*m);
    C1 = C*4*m;
    C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(1,:),y(1,:),n,l));
    C2 = C*(2*m+g*dt)*2*dt;
    x(2,n) = C1.*x(1,n) + C3(1) + C2.*v0(1,n);
    y(2,n) = C1.*y(1,n) + C3(2) + C2.*v0(2,n);
end

set(h,'XData',x(end,:),'YData',y(end,:));

drawnow;

%%

for t = 2:1000
    
    for n =1:N
        C = 1/(2*m-g*dt);
        C1 = C*4*m;
        C2 = -C*(2*m+g*dt);
        C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(t,:),y(t,:),n,l));
           
        x(t+1,n) = C1.*x(t,n) + C3(1) + C2*x(t-1,n);
        y(t+1,n) = C1.*y(t,n) + C3(2) + C2*y(t-1,n);
        
    end
    
    set(h,'XData',x(end,:),'YData',y(end,:));
    
%     for i = 1:N
%         addpoints(pl(i),t,x(end,i));
%     end
    
%     axis([t-20 20+t -0 10]);
    pause(0.1)
    drawnow;
    
    if KEY_IS_PRESSED
        close all;
        break;
    end
    
end


%%

function F = force(x,y,n,a)

F = [0 0];

if n-1 >= 1
    dx = x(n-1)-x(n);
    dy = y(n-1)-y(n);
    
    theta = atan(y/x);
    r = sqrt(dx^2+dy^2);
    delta = r-a;
    
    k = LenardJones(a,delta,1,2);
    
    F(1) = F(1) + k*(dx+delta*cos(theta));
    F(2) = F(2) + k*(dy+delta*sin(theta));
end

if n+1 <= length(x)
    dx = x(n)-x(n+1);
    dy = y(n)-y(n+1);
    
    theta = atan(y/x);
    r = sqrt(dx^2+dy^2);
    delta = r-a;
    
    k = LenardJones(a,delta,1,2);
    
    F(1) = F(1) - k*(dx+delta*cos(theta));
    F(2) = F(2) - k*(dy+delta*sin(theta));
end

end


function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
end
