clc;close all;clear all;
global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf;
set(gcf, 'KeyPressFcn', @myKeyPressFcn);
close all;

%% Pramètres de simulation

m = 1;
kb = 1.38e-23;
dt = 0.0001;
l = 1;

T0 = 10;
tmax = 1000000;

N=10;


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

global TYP
TYP = zeros(1,N);
v0 = zeros(2,N);
x = zeros(1,N);
y = zeros(1,N);

for n = 1:N
    TYP(n) = mod(n,2); 
    x(1,n) = randi(10);
    y(1,n) = randi(10);
    v0(1,n) = -1+2*rand(1);
    v0(2,n) = -1+2*rand(1);
end
% 
% x = [2 6 10  10 5];
% y = [2 10 9 6 2];
% v0 = [1 0.5 1 0.2 1; 0.3 0.3 -0.7 -0.6 0.5];


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
    
    T = temperature(x,y,N,kb,dt,m); % calcul de la température du système
    
    if t == 2
       T0 = T; 
       Tmax = T0*10;
    end
    
    %  Mise à jour des graphique à tout les 20 images
    if mod(t,100) == 0
        
        
        [xs,ys] = mc(x,y,N); % calcul du centre de masse
        
        addpoints(h2,t,T); % maj de T
        
        set(h1,'XData',x(end,:),'YData',y(end,:)); % maj des positions
        axis(f1,[-xx+xs xx+xs -yy+ys yy+ys]); % maj des axes
        
        drawnow;
    end
    
    t = t+1;
    
    %  Paramètre du thermostat
%     if T > tmax
%         
%         alp2 = sqrt(T0/T)^(-1);
%         v0 = zeros(2,N);
%         
%         for n =1:N
%             
%             vx = round((3*x(t,n) - 4*x(t-2,n)  + x(t-3,n))/(2*dt),3);
%             vy = round((3*y(t,n) - 4*y(t-2,n)  + y(t-3,n))/(2*dt),3);
% 
%             v0(:,n) = alp2 * [vx; vy];
%             
%         end
%         
%         for n = 1:N
%             
%             xi = normrnd(0,1);
%             
%             C = 1/(4*m);
%             C1 = C*4*m;
%             C3 = C*( 2*g*sqrt(2*D)*xi*dt^2 + 2*dt^2*force(x(1,:),y(1,:),n,l));
%             C2 = C*(2*m+g*dt)*2*dt;
%             x(t+1,n) = round(C1.*x(t,n) + C3(1) + C2.*v0(1,n),3);
%             y(t+1,n) = round(C1.*y(t,n) + C3(2) + C2.*v0(2,n),3);
%             
%         end
%         
%         t = t+1;
%         
%     end
    
    
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

global TYP

F = [0 0];

for i = 1:length(x)
    if i == n-1
        dx = x(i)-x(n);
        dy = y(i)-y(n);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        k = LenardJones(a,1,1,2);
        
        F(1) = F(1) + k*((delta)*cos(theta));
        F(2) = F(2) + k*((delta)*sin(theta));
        
    elseif i == n+1
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        
        k = LenardJones(a,1,1,2);
        
        F(1) = F(1) - k*((delta)*cos(theta));
        F(2) = F(2) - k*((delta)*sin(theta));
        
    end
    if i ~= n
        dx = x(n)-x(i);
        dy = y(n)-y(i);
        
        theta = atan2(dy,dx);
        
        r = sqrt(dx^2+dy^2);
        delta = r-a;
        
        if TYP(i) == TYP(n)
            e = 0.5;
        else
            e = -0.1;
        end
        
        f = LenardJones(a,delta,e,1);
        
        F(1) = F(1) + f*cos(theta);
        F(2) = F(2) + f*sin(theta);
        
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
