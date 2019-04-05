clc;clear all;close all;
%% Pramètres de simulation

dt = 1e-23; % pas de temps
T0 = 0.1; % Température à l'équilibre
N=2; % Longueur de la chaine
m = 10*1.66e-27; % masse d'un AA
a = 10e-10; % distance à l'équilibre entre les AAs

%THERMO = true;
THERMO_on = false;

tmax = 0.01e-18/dt;

fig_on = false;

%% Conditions initiale

global VOI ALT k

% Initialisation de variables
type = zeros(1,N);
v0 = zeros(2,N);
x = zeros(1000,N,'double');
y = zeros(1000,N,'double');

% Variable contenant la position des voisins d'un AA donné
VOI = diag(ones(1,N-1,'double'),1)+diag(ones(1,N-1,'double'),-1);
% Variable contenant la position des non-voisins d'un AA
ALT = ones(N,N)-VOI;
ALT = ALT - diag(ones(1,N));

for n = 1:N
    type(n) = mod(n,2); % type d'AA
    x(1,n) = 1.5*a*n; % position initiale en x
    y(1,n) = 0; % position initiale en y
    v = random_maxboltz(T0, 1, m); % grandeur de la vitesse
    ang = 360*rand(1); % orientation de la vitesse
    v0(1,n) = cosd(ang)*v; % vitesse initiale en x
    v0(2,n) = sind(ang)*v; % vitesse initiale en y
end

v0(1,:) = v0(1,:)-sum(v0(1,:));
v0(2,:) = v0(2,:)-sum(v0(2,:));


% Assigne la valeur de epsilon de l'interaction entre deux AAs
for i = 1:N
    for j = 1:N
        if type(i) == type(j)
            ALT(i,j) = ALT(i,j)* 1;
        else
            ALT(i,j) = ALT(i,j)* 0.2;
        end
    end
end

k =1e-4*LenardJones(a,1,1,2);

%% Test d'erreur

p = [2 4 8 16 32];

ex = zeros(length(p),2);
ey = zeros(length(p),2);
n = 1;

dt = 10e-24;
tmax = 0.1e-18/dt;

fprintf('Calcul avec dt = %.3e ... \n',dt)
tic
[x,y] = simulation(dt,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on);
temps = toc;
fprintf('Calcul exécuté en t = %.3f s \n',temps)

r = sqrt(x.^2+y.^2);

xmax = max(r-r(1,:));

for i = 1./p
    
    dtn = dt*i;
    tmax = 0.1e-18/dtn;
    
    fprintf('Calcul avec dt = 1/%d dt ... \n',p(n))
    tic
    [xn,yn] = simulation(dtn,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on);
    temps = toc;
    fprintf('Calcul exécuté en t = %.3f s \n',temps)
    
    rn = sqrt(xn.^2+yn.^2);
    
    ex(n,:) = abs(xmax - max(rn-rn(1,:)));
    %ex(n) = sum(abs(sum(r,1)./size(r,1) - sum(rn,1)./size(rn,1)));
    xmax = max(rn-rn(1,:));
    r=rn;
    
    n = n+1;

end

pft = polyfit(log(1./p),log(ex(:,1)'),1);
loglog(1./p,ex(:,1)');
