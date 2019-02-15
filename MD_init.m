clc; clear all; close all;
%% Param�tres 
global V epsi;

N       = 2;           % nombre d'AA
T       = 300;          % temp�rature
m       = 2.5e-25;      % masse d'un AA en Kg
g       = 0.5;            % gamma - coefficient de viscosit�
D       = 0.2;            % Coefficient de difusion
d       = 1;            % distance � potentiel null
xi     = 1;            % bruit thermique (bruit gaussien)
rayon   = 1;            % rayon

dt      = 0.001;         % pas de temps

duration = 100;         % temps de la simulation
fps      = 10;          % image par seconde
movie    = false;        % Cr�ation d'un fichier video si �gale true

%% Initialisation

r0      = 2^(1/6)*d;    % distance � l'�quilibre entre les AA li�
epsi = [1 0.2; 0.2 -1]; % matrice d'interaction entre les AA
type = randi(2,1,N); % type d'AA dans la chaine 2 pour hydrophile et 1 pour hydrophobe


V = [N T m g D d xi rayon r0 dt];

for n = 1:N
   protein(n) = AA(type(n), n,0);
end

MD_simulation(protein, duration, movie, fps)
