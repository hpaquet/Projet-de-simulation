
close all
rng('shufle')

%%%


A=[-1 0.5; 0.5 -1]; %% matrice des coefficients

T=0.5; % Temperature arbitraire
Hydrophile=2;
Hydrophobe=1;

%matrice des acides aminés
M=ones(8,8);
for i=1:8
    M(1,i)=2;
    M(i,1)=2;
    M(8,i)=2;
    M(i,8)=2;
end

%matrices des voisins cenvention [N,S,E,O]
P=cell(8);
P{1,1}=[0,2,2,0];
P{1,8}=[0,2,0,1];
P{8,1}=[2,0,2,0];
P{8,8}=[2,0,0,2];
for i=2:7
    P{1,i}=[0,1,2,2];
    P{i,1}=[2,2,1,0];
    P{8,i}=[1,0,2,2];
    P{i,8}=[2,2,0,1];
end
P{2,2}=[2,1,1,2];
P{2,7}=[2,1,2,1];
P{7,2}=[1,2,1,2];
P{7,7}=[1,2,2,1];
for i=3:6
    P{2,i}=[2,1,1,1];
    P{i,2}=[1,1,1,2];
    P{7,i}=[1,2,1,1];
    P{i,7}=[1,1,2,1];
end
for i=3:6
    for j=3:6
        P{i,j}=[1,1,1,1];
    end
end

%%%calcul Etot etqt initial
%Etot0=


%%boucle 
    %tirage aléatoire
    n=randi(8);
    m=randi(8);
            










