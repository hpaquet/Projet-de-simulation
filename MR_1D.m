clc;close all;clear all;
%%
% Fonctionne parfaitement
%
%
%
%%


global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);

%%

m = 1;

k = 10;

dt = 0.1;
l = 4;

N=2;

for n = 1:N
    x(1,n) = (n-randi([0,1]));
    v0(1,n) = 0;%(-2)^n;
end


subplot(2,1,1)
h=plot(0,0,'MarkerSize',100,'Marker','.');
axis([-8 15 -1 1]);
set(gca,'nextplot','replacechildren');

subplot(2,1,2)
for i =1:N
    pl(i) = animatedline('color',rand(1,3));
end

% pl(end+1) = animatedline('color',rand(1,3));
% pl(end+1) = animatedline('color',rand(1,3));

axis([0 100 -inf inf]);


%%
for n = 1:N
    x(2,n) = x(1,n) + dt^2*force(x(1,:),n,k,l)/(2*m) + dt.*v0(n) ;
end

set(h,'XData',x(end,:),'YData',zeros(1,N));

drawnow;
%%

for t = 2:1000
    
    for n =1:N
        x(t+1,n) = 2*x(t,n) - x(t-1,n) + dt^2*force(x(t,:),n,k,l)/m;    
    end
    
    set(h,'XData',x(end,:),'YData',zeros(1,N));
    
    for i = 1:N
        addpoints(pl(i),t,x(end,i));
    end
%     y1 = 1/2*(  cos(sqrt(2*k/m)*t*dt) + cos(-sqrt(2*k/m)*t*dt)   ) -1;
%     y2 = -1/2*(  cos(sqrt(2*k/m)*t*dt) + cos(-sqrt(2*k/m)*t*dt)  ) +3;
%     addpoints(pl(end-1),t, y1 );
%     addpoints(  pl(end),t, y2 );
    axis([t-20 20+t -inf inf]);
    %pause(0.1)
    drawnow;
    
    if KEY_IS_PRESSED
        close all;
        break;
    end
    
end




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


function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
    KEY_IS_PRESSED  = 1;
end
