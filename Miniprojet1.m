clc;close all;clear all;
%% Param�tres
R = 1; % rayon de la sph�re
n = 10; % nombre de dimension 
N_tot = 1e8; % nombre d'it�ration maximum 

Err0 = 1e-3; %erreur d�sir� 

%% Boucle principale
V_anal = (pi^(n/2))/(gamma(n/2 +1))*R^n; % Volume de l'hypershp�re obtenue de fa�on analytique 

range = 100:10:N_tot; 
err = [];
V = [];

% Calcul de volume de volume avec la m�thode de Mont�-Carlo 100 fois
for i = 1:100
    
    [V(end+1,:), err(end+1,:)] = volumeMC(R,n, V_anal, range);
    disp(i)
    
end

err_moy = sum(err,1)./100; % calcul de l'erreur moyenne


%% Question c) et d)

% calcul de l'erreur (pente du graphique loglog)
p = (log(err_moy(end))-log(err_moy(50)))/(log(range(end))-log(range(50)));

% valeur th�orique de p 
p_th = (log(1/sqrt(range(end)))-log(1/sqrt(range(50))))/(log(range(end))-log(range(50)));

% Valeur de A
A = err_moy./range.^p;
A = mean(A(round(length(A)/2):end));

% Calcul dla valeur de N pour avoir une erreur Err_0
N_perf = (Err0/A)^(1/p);

%% Graphique

range2 = 100:10:N_perf; 

loglog(range,err_moy')
hold on
loglog(range2, A*range2.^p,'--r')
loglog(N_perf,A*N_perf^p ,'ok')
hold off
xlabel('N')
ylabel('moyenne |V(N)-V_anal|/V_anal')

%Affichage d'un message avec les r�sultats
m0 = sprintf('\bfbold\rm Hypersph�re de dimension %2d \n', n );
m1 = sprintf('L''erreur (p) est de : %.3f \n', p );
m2 = sprintf('La valeur th�orique de p est de : %.3f \n', p_th );
m3 = sprintf('La valeur de N qui garantit (en moyenne) une precision relative Err_0  est de : %.3f \n', N_perf );
msgbox({m0, m1, m2, m3},'R�sultat du calcul','warn')


%% Calcul du volume � l'aide de la m�thode de Mont� Carlo

function [V, err] = volumeMC(R,n, Vr, range)

    V = [];
    err = [];
    N_int = 0;
    N_tot = 0;
    
    for k = range
        
        point = 2*R*rand(1,n)-R; % Cr�ation d'un point (vecteur de coordonn�e) dans la boite de cot� 2R
        N_tot = N_tot +1;
        
        if (norm(point))^2 <= R^2 % v�rification si le point est dans l'hyp�rshpere
            N_int = N_int+1;
        end

        V(end+1) = N_int/N_tot * 2^n; % Calcul du volume 
        err(end+1) = norm((V(end)-Vr)/Vr); % Calcul de l'erreur
        
    end
    
end