clc; clear all; close all;
%% Param�tres 
global V epsi;

N       = 2;           % nombre d'AA
T       = 300;          % temp�rature
m       = 2.5e-25;      % masse d'un AA en Kg
g       = 1;            % gamma - coefficient de viscosit�
D       = 0.2;            % Coefficient de difusion
a       = 1;            % distance � potentiel null
xi     = 0;            % bruit thermique (bruit gaussien)
rayon   = 1;            % rayon

dt      = 0.01;         % pas de temps

duration = 1000;         % temps de la simulation
fps      = 10;          % image par seconde
movie    = false;        % Cr�ation d'un fichier video si �gale true

%% Initialisation

epsi = [1 0.2; 0.2 -1]; % matrice d'interaction entre les AA
type = randi(2,1,N); % type d'AA dans la chaine 2 pour hydrophile et 1 pour hydrophobe


V = [N T m g D a xi rayon dt];

protein = zeros(duration,2*N);

for n = 1:N
   protein(1,2*n-1) = n*a;
end

%v0 = random_maxboltz(T,N,m);

v0 = zeros(1,N);


MD_simulation(protein, type, duration, v0);
