clc;close all;clear all;
%%

global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);

%%

m = 1;
k = 1;

dt = 0.1;
l = 2;

N=2;

for n = 1:N
    x(1,n) = n;
    y(1,n) = n;
    v0(1,n) = 0;%(-1)^(n-1);%(-1)^n;
    v0(2,n) = 0;
end


% subplot(2,1,1)

h=plot(0,0,'MarkerSize',100,'Marker','.','LineWidth',5);
axis([-10 20 -5 5]);
set(gca,'nextplot','replacechildren');

% subplot(2,1,2)
% for i =1:N
%     pl(i) = animatedline('color',rand(1,3));
% end
% axis([0 100 -0 10]);
%%

for n = 1:N
    
    C3 = dt^2*force(x(1,:),y(1,:),n,k,l)/(2*m);

    x(2,n) = x(1,n) + C3(1) + dt*v0(1,n);
    y(2,n) = y(1,n) + C3(2) + dt*v0(2,n);
    
end

set(h,'XData',x(end,:),'YData',y(end,:));

drawnow;

%%

for t = 2:10000
    
    for n =1:N
        
        C3 = dt^2*force(x(t,:),y(t,:),n,k,l)/m;
        
        x(t+1,n) = 2*x(t,n) + C3(1) - x(t-1,n);
        y(t+1,n) = 2*y(t,n) + C3(2) - y(t-1,n);
        
    end
    
    set(h,'XData',x(end,:),'YData',y(end,:));
    
    %     for i = 1:N
    %         addpoints(pl(i),t,x(end,i));
    %     end
    
    %     axis([t-20 20+t -0 10]);
    %pause(0.1)
    drawnow;
    
    if KEY_IS_PRESSED
        close all;
        break;
    end
    
end


%%

function F = force(x,y,n,k,a)

F = [0 0];

if n-1 >= 1
    
    dx = x(n)-x(n-1);
    dy = y(n)-y(n-1);
    
    theta = atan2(dy,dx);
    
    r = sqrt(dx^2+dy^2);
    delta = r-a;
    
    F(1) = F(1) - k*((delta)*cos(theta));
    F(2) = F(2) - k*((delta)*sin(theta));
    
end

if n+1 <= length(x)
    
    dx = x(n)-x(n+1);
    dy = y(n)-y(n+1);
    
    theta = atan2(dy,dx);
    
    r = sqrt(dx^2+dy^2);
    delta = r-a;
    
    F(1) = F(1) - k*((delta)*cos(theta));
    F(2) = F(2) - k*((delta)*sin(theta));
    
end

end


function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
end
