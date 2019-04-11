clc;clear all;close all;
%% Pramètres de simulation


T0 = 1e20; % Température à l'équilibre
N=64; % Longueur de la chaine
m = 10*1.66e-27; % masse d'un AA
a = 10e-10; % distance à l'équilibre entre les AAs

%THERMO = true;
THERMO_on = false;

%%

fig_on = true; % active l'affichage

test = false; % test de convergence

simu = true; % simulation infinie avec affichage

repl = false; % test du repliement

testre = false;

%% Conditions initiale

global VOI ALT k type

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
    type(n) = 1; % type d'AA
    x(1,n) = 1.5*a*n; % position initiale en x
    y(1,n) = 0; % position initiale en y
    v = random_maxboltz(T0, 1, m); % grandeur de la vitesse
    ang = 360*rand(1); % orientation de la vitesse
    v0(1,n) = cosd(ang)*v; % vitesse initiale en x
    v0(2,n) = sind(ang)*v; % vitesse initiale en y
end

type = [1 1 1 1 0 0 1 1 0 0 1 1 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1];

% type = [ones(1,4) zeros(1,2) ones(1,2) zeros(1,22) ones(1,20)  zeros(1,22)];
v0(1,:) = v0(1,:)-sum(v0(1,:))/N;
v0(2,:) = v0(2,:)-sum(v0(2,:))/N;


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

k =100e-4*LenardJones(a,1,1,2);


%% Simulation normal
if simu == true
    dt = 10e-24; % pas de temps
   [xn,yn] = simulation(dt,0.01/dt,T0,N,m,a,false,v0,x,y,true); 
end

%% Test d'erreur

if test == true    
    
p = [2 4 8 16];

ex = zeros(length(p),N);
ey = zeros(length(p),2);
n = 1;

dt = 10e-24;
tmax = 0.01e-18/dt;

fprintf('Calcul avec dt = %.3e ... \n',dt)
tic
[x,y] = simulation(dt,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on);
temps = toc;
fprintf('Calcul exécuté en t = %.3f s \n',temps)

r = sqrt(x.^2+y.^2);

for i = p
    
    dtn = dt/i;
    tmax = 0.01e-18/dtn;
    
    fprintf('Calcul avec dt = 1/%d dt ... \n',p(n))
    tic
    [xn,yn] = simulation(dtn,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on);
    temps = toc;
    fprintf('Calcul exécuté en t = %.3f s \n',temps)
    
    rn = sqrt(xn.^2+yn.^2);
    
    ex(n,:) = max(abs(r(1:end,:) - rn(1:2:end,:)));

    r=rn;
    
    n = n+1;

end

pft = polyfit(log(1./p),log(ex(:,1)'),1);
loglog(1./p,ex(:,1)');

end

%% Test d'erreur

if testre == true    
    
p = [2 4 8 16 32];

ex = zeros(length(p),N);
ey = zeros(length(p),N);

n = 1;

dt = 10e-24;
t = 0.01e-18;
tmax = t/dt;


for i = p
    
    dtn = dt/i;
    tmax = t/dtn;
    
    fprintf('Calcul avec dt = 1/%d dt ... \n',p(n))
    tic
    [xn,yn] = simulation(dtn,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on);
    temps = toc;
    fprintf('Calcul exécuté en t = %.3f s \n',temps)
    
    rn = sqrt(xn.^2+yn.^2);
    
    xr = x(1,:)'*cos(sqrt(k/m)*(dtn:dtn:t))+v0(1,:)'/sqrt(k/m)*sin(sqrt(k/m)*(dtn:dtn:t));
    yr = v0(2,:)'/sqrt(k/m)*sin(sqrt(k/m)*(dtn:dtn:t));
    r = sqrt(xr.^2+yr.^2);
    
    ex(n,:) = sum(abs(xr' - xn),1);
    ey(n,:) = sum(abs(yr' - yn),1);
    

    %r=rn;
    
    n = n+1;

end

pft = polyfit(log(1./p),log(ex(:,1)'),1);
loglog(1./p,ex(:,1)');

end



%% Test de repliement

if repl == true
    
    dt = 100e-24;
    tmax = 1e-18/dt;
    
    xmax = [];xmin = [];xs = [];ys = [];ymax = [];ymin = [];
    
    l = [];h = [];

    figure()
    hold on
    
    for i = 1:5
        
        for n = 1:N
            x(1,n) = 1.5*a*n; % position initiale en x
            v = random_maxboltz(T0, 1, m); % grandeur de la vitesse
            ang = 360*rand(1); % orientation de la vitesse
            v0(1,n) = cosd(ang)*v; % vitesse initiale en x
            v0(2,n) = sind(ang)*v; % vitesse initiale en y
        end
        
        v0(1,:) = v0(1,:)-sum(v0(1,:))/N;
        v0(2,:) = v0(2,:)-sum(v0(2,:))/N;
        
        [x,y] = simulation(dt,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on);
       
        % Centre de masse de la protéine
        xs(end+1) = sum(x(end,:))/N;
        ys(end+1) = sum(y(end,:))/N;
        
        % interval des AAs
        l = max(x(end,:)) - min(x(end,:));
        h = max(y(end,:)) - min(y(end,:));
        
        a = -h:1e-11:h;
        b = l/h.*sqrt(-x.^2+h^2);
        
        plot(a,b)
        
    end

    hold off
    
    
end





