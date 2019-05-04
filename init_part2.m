clc;clear all;close all;
%% Pramètres de simulation


T0 = 1e20; % Température à l'équilibre
N=64; % Longueur de la chaine
m = 75*1.66e-27; % masse d'un AA
a = 10e-10; % distance à l'équilibre entre les AAs

%THERMO = true;
THERMO_on = false;

%%

fig_on = false; % active l'affichage

test = false; % test de convergence

simu = true; % simulation infinie avec affichage

repl = false; % test du repliement

testre = false;

stab  = false; % Analyse de stabilité

movie = false;

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
    type(n) = mod(n,2); % type d'AA
    x(1,n) = 1*a*n; % position initiale en x
    y(1,n) = 0; % position initiale en y
    v = random_maxboltz(T0, 1, m); % grandeur de la vitesse
    ang = 360*rand(1); % orientation de la vitesse
    v0(1,n) = cosd(ang)*v; % vitesse initiale en x
    v0(2,n) = sind(ang)*v; % vitesse initiale en y
end

% type = [ones(1,4) zeros(1,2) ones(1,2) zeros(1,22) ones(1,20)  zeros(1,22)];
v0(1,:) = v0(1,:)-sum(v0(1,:))/N;
v0(2,:) = v0(2,:)-sum(v0(2,:))/N;

% Assigne la valeur de epsilon de l'interaction entre deux AAs
for i = 1:N
    for j = 1:N
        if (type(i) == type(j) ) && type(i) == 1 
            ALT(i,j) = -ALT(i,j)*0.99;
        elseif (type(i) == type(j) ) && type(i) == 0
            ALT(i,j) = ALT(i,j)*0.99;
        else
            ALT(i,j) = -ALT(i,j)* 0.5;
        end
    end
end

k =1e1*LenardJones(a,1,1,2);


%% Simulation normal
if simu == true
    dt = 3e-22; % pas de temps
        type = [1 1 1 1 0 0 1 1 ...
            0 0 1 1 0 0 1 1 ...
            1 1 1 0 0 0 0 0 ...
            0 0 0 0 0 0 0 1 ...
            1 0 0 0 0 0 0 1 ...
            1 0 0 0 0 0 0 1 ...
            1 0 0 0 0 0 0 1 ...
            1 1 1 1 1 1 1 1];
    [xn,yn] = simulation(dt,1e-18/dt,T0,N,m,a,THERMO_on,v0,x,y,true,movie);
end

%% Test d'erreur

if test == true
    
    p = [2 4 8 16 32 64 128];
    
    ex = zeros(length(p),N);
    ey = zeros(length(p),2);
    n = 1;
    
    dt = 10e-24;
    tmax = 0.01e-18/dt;
    
    fprintf('Calcul avec dt = %.3e ... \n',dt)
    tic
    [x,y] = simulation(dt,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on,false);
    temps = toc;
    fprintf('Calcul exécuté en t = %.3f s \n',temps)
    
    r = sqrt(x.^2+y.^2);
    
    for i = p
        
        dtn = dt/i;
        tmax = 0.01e-18/dtn;
        
        fprintf('Calcul avec dt = 1/%d dt ... \n',p(n))
        tic
        [xn,yn] = simulation(dtn,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on,false);
        temps = toc;
        fprintf('Calcul exécuté en t = %.3f s \n',temps)
        
        rn = sqrt(xn.^2+yn.^2);
        
        ex(n,:) = max(abs(r(1:end,:) - rn(1:2:end,:)));
        
        r=rn;
        
        n = n+1;
        
    end
    
    pft = polyfit(log(1./p),log(ex(:,1)'),1);
    loglog(dt./p,ex(:,1)');
    title('Erreur en fonction du pas de temps')
    text(1.8e-26,0.3e-14,'Ordre de convergence (pente) = 1.937')
    xlabel('Log(dt)')
    ylabel('Log(Erreur)')
    
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
        [xn,yn] = simulation(dtn,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on,false);
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
        
        [x,y] = simulation(dt,tmax,T0,N,m,a,THERMO_on,v0,x,y,fig_on,false);
        
        % Centre de masse de la protéine
        xs(end+1) = sum(x(end,:))/N;
        ys(end+1) = sum(y(end,:))/N;
        
        % interval des AAs
        l = max(x(end,:)) - min(x(end,:));
        h = max(y(end,:)) - min(y(end,:));
        

    end
    
    hold off
    
    
end

%% Analyse de stabilité

if stab == true
    
    a = 10e-10;
    k =1e1*LenardJones(a,1,1,2);
    
    type = [1 1 1 1 0 0 1 1 ...
            0 0 1 1 0 0 1 1 ...
            1 1 1 0 0 0 0 0 ...
            0 0 0 0 0 0 0 1 ...
            1 0 0 0 0 0 0 1 ...
            1 0 0 0 0 0 0 1 ...
            1 0 0 0 0 0 0 1 ...
            1 1 1 1 1 1 1 1];
    
    x(1,:) = a*[ 3 2 1 1 2 2 1 1 ...
                2 2 1 1 2 2 1 1 ...
                2 3 4 4 3 3 4 4 ...
                3 3 4 4 3 3 4 4 ...
                5 5 5 5 5 5 5 5 ...
                6 6 6 6 6 6 6 6 ...
                7 7 7 7 7 7 7 7 ...
                8 8 8 8 8 8 8 8 ];
    y(1,:) = a*[ 1 1 1 2 2 3 3 4 ...
                4 5 5 6 6 7 7 8 ...
                8 8 8 7 7 6 6 5 ...
                5 4 4 3 3 2 2 1 ...
                1 2 3 4 5 6 7 8 ...
                8 7 6 5 4 3 2 1 ...
                1 2 3 4 5 6 7 8 ...
                8 7 6 5 4 3 2 1 ];
            
    dt = 1e-23; % pas de temps
    err = [];
    
    for i = 1:20
        T0n = T0*10^i;
        for n = 1:N
            v = random_maxboltz(T0n, 1, m); % grandeur de la vitesse
            ang = 360*rand(1); % orientation de la vitesse
            v0(1,n) = cosd(ang)*v; % vitesse initiale en x
            v0(2,n) = sind(ang)*v; % vitesse initiale en y
        end
        v0(1,:) = v0(1,:)-sum(v0(1,:))/N;
        v0(2,:) = v0(2,:)-sum(v0(2,:))/N;
        
        [xn,yn] = simulation(dt,0.02e-18/dt,T0n,N,m,a,false,v0,x,y,fig_on,false);
        
        rn = sqrt(xn.^2+yn.^2);
        
        err(end+1) = sum(abs(rn(1,:)-rn(end,:)));
        
        disp(i)
        
    end
    
    loglog(T0.*10.^(1:20),err)
    title('Erreur en fonction de la température')
    xlabel('Log de la température')
    ylabel('Log de l''erreur')
    
    
end




