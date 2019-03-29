clc;close all;clear all;

%% Pramètres de simulation

m = 100*1.66e-27; % masse d'un AA
a = 5e-8; % distance à l'équilibre entre les AAs

dt = 1e-25; % pas de temps
T0 = 10; % Température à l'équilibre

N=10; % Longueur de la chaine


%% Paramètre graphique

% Dimension du graphique 
xx = 4e-8; 
yy = 4e-8;

figure('Name','Dynamique moléculaire','Position', [ 50 50 1200 600 ],'NumberTitle','off');
f1 = subplot(1,2,1);
h1=plot(0,0,'MarkerSize',50,'Marker','.','LineWidth',3);
axis off
axis(f1,[-xx xx -yy yy])
set(gca,'nextplot','replacechildren')
subplot(1,2,2)
h2=animatedline;
ylabel('Température (K)')
xlabel('Temps (as)')

%% Conditions initiale

global VOI ALT

% Initialisation de variables
type = zeros(1,N);
v0 = zeros(2,N);
x = zeros(1000,N,'single');
y = zeros(1000,N,'single');

% Variable contenant la position des voisins d'un AA donné
VOI = diag(ones(1,N-1,'single'),1)+diag(ones(1,N-1,'single'),-1);
% Variable contenant la position des non-voisins d'un AA
ALT = zeros(N,N)-VOI;


for n = 1:N
    type(n) = mod(n,2); % type d'AA 
    x(1,n) = rand(1)*20e-9; % position initiale en x
    y(1,n) = rand(1)*20e-9; % position initiale en y
    v = random_maxboltz(T0, 1, m); % grandeur de la vitesse 
    ang = pi*rand(1); % orientation de la vitesse
    v0(1,n) = cos(ang)*v; % vitesse initiale en x
    v0(2,n) = sin(ang)*v; % vitesse initiale en y
end

% Assigne la valeur de epsilon de l'interaction entre deux AAs
for i = 1:N
    for j = 1:N
        if type(i) == type(n)
            ALT(i,j) = ALT(i,j)* 1;
        else
            ALT(i,j) = ALT(i,j)* -0.5;
        end
    end
end


%% Première itération (t = 1)

C3 = dt^2*force(x(1,:),y(1,:),a,N)/(2*m);

x(2,:) = x(1,:) + C3(1,:) + dt*v0(1,:); 
y(2,:) = y(1,:) + C3(2,:) + dt*v0(2,:);
    

%% Boucle principale (pour t>1)

t = 2;

while(true)
    
    % Mise à jour des positions
        
    C3 = dt^2*force(x(t,:),y(t,:),a,N)/m;

    x(t+1,:) = 2*x(t,:) + C3(1,:) - x(t-1,:);
    y(t+1,:) = 2*y(t,:) + C3(2,:) - y(t-1,:);
        
    t = t+1;
    
    % Température du système
    T = temperature(x,y,N,dt,m,t);
    
    %  Mise à jour des graphique à tout les 10 images
    if mod(t,50) == 0
        
        [xs,ys] = mc(x,y,N,t); % calcul du centre de masse
        
        addpoints(h2,t*dt/1e-18,T); % maj du graphe de la temperature
        
        set(h1,'XData',x(t,:),'YData',y(t,:)); % maj des positions
        axis(f1,[-xx+xs xx+xs -yy+ys yy+ys]); % maj des axes
        
        drawnow;
    end
    
    %  Thermostat 
    if (T > 10*T0 || T< T0/10 ) && t>3
        
        alp2 = sqrt(T0/T); % coefficient alpha 
            
        % Valeur des vitesses
        vx = (3*x(t,:) - 4*x(t-2,:)  + x(t-3,:))./(2*dt);
        vy = (3*y(t,:) - 4*y(t-2,:)  + y(t-3,:))./(2*dt);

        % Ajustement des vitesses selon alpha
        v0 = alp2 * [vx; vy];
        
        % Première itération de la simulation pour repartir la simulation
        C3 = dt^2*force(x(1,:),y(1,:),a,N)/(2*m);
        
        x(t+1,:) = x(t,:) + C3(1,:) + dt*v0(1,:);
        y(t+1,:) = y(t,:) + C3(2,:) + dt*v0(2,:);
        
        t = t+1;
        
    end
    
    % initialise 1000 nouveau espace mémoire dans les variable de position
    if mod(t,1000) == 0
        x = [x;zeros(1000,N)];
        y = [y;zeros(1000,N)];
    end
    
    profile on
    
end

profile viewer

%% Fonctions

% Calcul le centre de masse de la protéine
function [xs,ys] = mc(x,y,N,t)

xs = sum(x(t,:))/N;
ys = sum(y(t,:))/N;

end


% Calcul la force entre les AAs
function F = force(x,y,a,N)

global VOI ALT

di = eye(N)>0;
k = LenardJones(a,1,20,2); % Constante de rapel

dx = (x.*ones(N,N))-x';dy = (y.*ones(N,N))-y'; % Distance x et y entre les AAs
r = sqrt(dx.^2+dy.^2); % Distance entre les AAs
delta = r-a; % Distance du point d'équilibre des AAs

% Force entre les voisins (liaisons covalente)
Fx = VOI.*k.*delta.*dx./r;Fy = VOI.*k.*delta.*dy./r;
Fx(di) = 0;Fy(di) = 0;
Fx = sum(Fx,2);Fy = sum(Fy,2);
F = [Fx';Fy'];

% Force entre les autres
Fx = ALT.*LenardJones(a,delta,ALT,1).*dx./r;Fy = ALT.*LenardJones(a,delta,ALT,1).*dx./r;
Fx(di) = 0;Fy(di) = 0;
Fx = sum(Fx,2);Fy = sum(Fy,2);
F = F+[Fx';Fy'];

end

% Calcule la température du système
function T = temperature(x,y,N,dt,m,t)

C = m/(8*1*N*dt^2);
%C = m/(8*physconst('Boltzmann')*N*dt^2);

vx = 3*x(t,:)-4*x(t-1,:)+x(t-2,:);
vy = 3*y(t,:)-4*y(t-1,:)+y(t-2,:);
v2 = vx.^2+vy.^2;

T = double(C*sum(v2));


end
