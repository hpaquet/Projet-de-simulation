clc;close all; clear all;
%% Question iii)
a = 0.16/100; % alpha
% a = 0.1/100; % alpha
m = 50/1000; % masse
k = 20/100; % constante de ressort

% calculs préliminaire
T = 2*pi*sqrt(m/k);
t_max = 0.6*T;

expo = 3:0.25:4.75;
delt = 10.^(-expo)*T/pi;

% initialisation de variable
d_max = [];
Err = [];
j = 1;

% Calculs 

for dt = delt
    for i = [1 2]
    
    [y,time] = simu(dt/i,a,t_max,k,m);
    
    d = abs(y-y(1));
    d_max(j,i) = max(d);
    
    end
    
    Err(j) = abs(d_max(j,1)-d_max(j,2));
    j = j+1;
    
end

% Affichage
figure()
subplot()
plot(delt,Err)
title('Erreur en fonction du pas de temps')
xlabel('Pas de temps (dt)')
ylabel('Erreur')


%% Question iiii)

dt = 1e-3 *T/pi;
t_max = 50*T;
yeq = -2.45;
[y,time] = simu(dt,a,t_max,k,m);

y = y - yeq;

p = findpeaks(y,'MinPeakHeight',yeq);
p = length(p);

x = ones(1,length(time))*max(y)/10;
figure()
hold on
plot(time,y(1:end-1))
plot(time,x)
hold off
title('Position par rapport au point d''équilibre en fonction du temps')
xlabel('Temps (t)')
ylabel('Position (cm)')

%% Résolution de l'équation différentielle 

function [y, time] = simu(dt,a,t_max,k,m)

    % paramètres de résolution
    g = 9.8;
    v0 = -2.5;
    y = -0.5;
    
%     v0 = -1.5;
%     y = -1;

    % initialisation des valeurs
    
    time = 0:dt:t_max;
    
    C1 = (2+a*dt/m)^(-1);
    C2 = 4-2*k*dt^2/m;
    C3 = a*dt/m-2;
    C4 = 2*g*dt^2;
    
    % Résolution
    
    for t = time+dt

        if t == dt % première itération
            y(end+1) = y(1)*( 1 - k/(2*m)*dt^2 ) + v0*dt*( 1-a/(2*m)*dt ) - dt^2*g/2;
        else
            y(end+1)= C1*( C2*y(end) + C3*y(end-1) - C4 );
        end

    end
end