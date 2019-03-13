clc;close all;clear all;
global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);
close all;

%% Pramètres de simulation

m = 1;

k = 1;

kb = 1.38e-23;
dt = 0.0006;
l = 1;
g = 1e-20;
xi = 1;
D = 1/g;
alp2 = 1;

T0 = 1;
tmax = 2;

N=3;


%% Paramètre graphique

T = [];
xx = 10;
yy = 10;

figure('Name','Dynamique moléculaire','Position', [ 50 50 1200 600 ],'NumberTitle','off');
f1 = subplot(1,2,1);
h1=plot(0,0,'MarkerSize',100,'Marker','.','LineWidth',5);
axis off
axis(f1,[-xx xx -yy yy])
set(gca,'nextplot','replacechildren')
subplot(1,2,2)
h2=animatedline;


%% Conditions initiale

for n = 1:N
    x(1,n) = randi(10);
    y(1,n) = randi(10);
    v0(1,n) = -1+2*rand(1);
    v0(2,n) = -1+2*rand(1);
end

x = [9 7 4 ];
y = [8 9 6];
v0 = [-0.6 0.5 -0.5; -0.2 1 0.7];


%% Première itération (t = 1)

for n = 1:N
    
    xi = normrnd(0,1);
    
    C = 1/(4*m);
    C1 = C*4*m;
    C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(1,:),y(1,:),n,l));
    C2 = C*(2*m+g*dt)*2*dt;
    x(2,n) = round(C1.*x(1,n) + C3(1) + C2.*v0(1,n),3);
    y(2,n) = round(C1.*y(1,n) + C3(2) + C2.*v0(2,n),3);
    
end


%% Boucle principale (pour t>1)

t = 2;

while(true)
    
    % Mise à jour des positions
    for n =1:N
        
        xi = normrnd(0,1);
        
        C = 1/(2*m-alp2*g*dt);
        C1 = C*4*m;
        C2 = -C*(2*m+alp2*g*dt);
        C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(t,:),y(t,:),n,l));
        
        x(t+1,n) = round(C1.*x(t,n) + C3(1) + C2*x(t-1,n),3);
        y(t+1,n) = round(C1.*y(t,n) + C3(2) + C2*y(t-1,n),3);
        
    end
    
    
    %  Mise à jour des graphique à tout les 20 images
    if mod(t,100) == 0
        
        T = temperature(x,y,N,kb,dt,m); % calcul de la température du système
        [xs,ys] = mc(x,y,N); % calcul du centre de masse
        
        addpoints(h2,t,T); % maj de T
        
        set(h1,'XData',x(end,:),'YData',y(end,:)); % maj des positions
        axis(f1,[-xx+xs xx+xs -yy+ys yy+ys]); % maj des axes
        
        drawnow;
        
    end
    
    
    % Paramètre du thermostat
    if T > tmax
        alp2 = sqrt(T0/T);
        for n =1:N
         
            x(t+1,n) = round(alp2*(x(t,n) - x(t-2,n) ) + x(t-1,n),3);
            y(t+1,n) = round(alp2*(y(t,n) - y(t-2,n) ) + y(t-1,n),3);
            
        end
        alp2 = 1;
    else; alp2 = 1;
    end
    
    t = t+1;
    
    
    
    if KEY_IS_PRESSED
        close all;
        break;
    end
end


%% Fonctions

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
        
        f = LenardJones(a,delta,-10,1);
        
        F(1) = F(1) - f*cos(theta);
        F(2) = F(2) - f*sin(theta);
        
    end
end
end


function T = temperature(x,y,N,k,dt,m)

T = 0;
C = m/(8*1*N*dt^2);
%C = m/(8*k*N*dt^2);

for i = 1:N
    
    vx = 3*x(end,i)-4*x(end-1,i)+x(end-2,i);
    vy = 3*y(end,i)-4*y(end-1,i)+y(end-2,i);
    
    v2 = vx^2+vy^2;
    
    T = T + v2;
    
end

T = C*T;

end

%%

function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
end
