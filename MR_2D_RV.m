clc;close all;clear all;
global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);
close all;

%% Pramètres de simulation

m = 100*1.66e-27;
kb = 1.38e-23;
dt = 1e-25;
l = 10e-8;

T0 = 0.5;
tmax = 30;

N=2;


%% Paramètre graphique

T = [];
xx = 10e-9;
yy = 10e-9;

figure('Name','Dynamique moléculaire','Position', [ 50 50 1200 600 ],'NumberTitle','off');
f1 = subplot(1,2,1);
h1=plot(0,0,'MarkerSize',50,'Marker','.','LineWidth',3);
axis off
axis(f1,[-xx xx -yy yy])
set(gca,'nextplot','replacechildren')
subplot(1,2,2)
h2=animatedline;
ylabel('Température (T)')
xlabel('Temps (t)')

%% Conditions initiale

global TYP
TYP = zeros(1,N);
v0 = zeros(2,N);
x = zeros(1000,N);
y = zeros(1000,N);

for n = 1:N
    TYP(n) = mod(n,2);
    x(1,n) = rand(1)*10e-9;
    y(1,n) = rand(1)*10e-9;
    v0(1,n) = (-1+2*rand(1))*random_maxboltz(T0, 1, m);
    v0(2,n) = (-1+2*rand(1))*random_maxboltz(T0, 1, m);
end


%% Première itération (t = 1)


for n = 1:N
    
    C3 = dt^2*force(x(1,:),y(1,:),n,l)/(2*m);
    
    x(2,n) = x(1,n) + C3(1) + dt*v0(1,n);
    y(2,n) = y(1,n) + C3(2) + dt*v0(2,n);
    
end


%% Boucle principale (pour t>1)

t = 2;

while(true)
    
    
    % Mise à jour des positions
    for n =1:N
        
        C3 = dt^2*force(x(t,:),y(t,:),n,l)/m;
        
        x(t+1,n) = 2*x(t,n) + C3(1) - x(t-1,n);
        y(t+1,n) = 2*y(t,n) + C3(2) - y(t-1,n);
        
    end
    
    T = temperature(x,y,N,kb,dt,m,t+1); % calcul de la température du système
    
    %  Mise à jour des graphique à tout les 20 images
    if mod(t,100) == 0
        
        [xs,ys] = mc(x,y,N,t+1); % calcul du centre de masse
        
        addpoints(h2,t,T); % maj de T
        
        set(h1,'XData',x(t+1,:),'YData',y(t+1,:)); % maj des positions
        axis(f1,[-xx+xs xx+xs -yy+ys yy+ys]); % maj des axes
        
        drawnow;
    end
    
    t = t+1;
    
    %  Paramètre du thermostat
    if (T > 10*T0 || T< T0/10 ) && t>3
        
        alp2 = sqrt(T0/T);
        v0 = zeros(2,N);
        
        for n =1:N
            
            vx = (3*x(t,n) - 4*x(t-2,n)  + x(t-3,n))/(2*dt);
            vy = (3*y(t,n) - 4*y(t-2,n)  + y(t-3,n))/(2*dt);
            
            v0(:,n) = alp2 * [vx; vy];
            
        end
        
        for n = 1:N
            
            C3 = dt^2*force(x(t,:),y(t,:),n,l)/(2*m);
            
            x(t+1,n) = x(t,n) + C3(1) + dt*v0(1,n);
            y(t+1,n) = y(t,n) + C3(2) + dt*v0(2,n);
            
        end
        
        t = t+1;
        
    end
    
    if mod(t,1000) == 0
        x = [x;zeros(1000,N)];
        y = [y;zeros(1000,N)];
    end
    
    
    if KEY_IS_PRESSED
        close all;
        break;
    end
    
    
end


%% Fonctions

function [xs,ys] = mc(x,y,N,t)
xs = sum(x(t,:))/N;
ys = sum(y(t,:))/N;
end


function F = force(x,y,n,a)

global TYP

F = [0 0];

for i = 1:length(x)
    if i == n-1
        dx = x(i)-x(n);
        dy = y(i)-y(n);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        k = LenardJones(a,1,20,2);
        
        F(1) = F(1) + k*((delta)*(dx/r));
        F(2) = F(2) + k*((delta)*(dy/r));
        
    elseif i == n+1
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        k = LenardJones(a,1,20,2);
        
        F(1) = F(1) - k*((delta)*(dx/r));
        F(2) = F(2) - k*((delta)*(dy/r));
        
    elseif i ~= n
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        if TYP(i) == TYP(n)
            e = 1;
        else
            e = -0.5;
        end
        
        f = LenardJones(a,delta,e,1);
        
        F(1) = F(1) + f*dx/r;
        F(2) = F(2) + f*dy/r;
        
    end
end
end


function T = temperature(x,y,N,k,dt,m,t)

T = 0;
C = m/(8*1*N*dt^2);
%C = m/(8*k*N*dt^2);

for i = 1:N
    
    vx = 3*x(t,i)-4*x(t-1,i)+x(t-2,i);
    vy = 3*y(t,i)-4*y(t-1,i)+y(t-2,i);
    
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
