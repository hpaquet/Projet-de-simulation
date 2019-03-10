
close all
rng('shufle')

%%%


A=[-1 0.5; 0.5 -1]; %% matrice des coefficients

T=0.5; % Temperature arbitraire
Hydrophile=2;
Hydrophobe=1;

%matrice des acides aminés
M_AA=ones(8,8);
for i=1:8
    M_AA(1,i)=2;
    M_AA(i,1)=2;
    M_AA(8,i)=2;
    M_AA(i,8)=2;
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
%Etot0 = 1


%%boucle 
    %tirage aléatoire
    n=randi(8);
    m=randi(8);          
%%changement d'AA  et changement des voisins pour ses voisins
  %%changement de L'AA
  if (M(m,n)==1)
    M(m,n)=2;
  else
    M(m,n)=2;
  end
  %changements pour les voisins
  if (m==1)
    if(n==1)
      P{m+1,n}{1}= M(m,n); %voisin inférieur
      P{m,n+1}{4}= M(m,n); %voisin de droite
    else if(n==8)
      P{m+1,n}{1}= M(m,n); %voisin inférieur
      P{m,n-1}{3}= M(m,n); %voisin de gauche
    else
      P{m+1,n}{1}= M(m,n); %voisin inférieur
      P{m,n-1}{3}= M(m,n); %voisin de gauche
      P{m,n+1}{4}= M(m,n); %voisin de droite
    end
  else if (m==8)
     if (n==1)
       P{m-1,n}{2}= M(m,n); %voisin supérieur
       P{m,n+1}{4}= M(m,n); %voisin de droite
     else if (n==8)
       P{m-1,n}{2}= M(m,n); %voisin supérieur
       P{m,n-1}{3}= M(m,n); %voisin de gauche
     else
       P{m-1,n}{2}= M(m,n); %voisin supérieur
      P{m,n-1}{3}= M(m,n); %voisin de gauche
      P{m,n+1}{4}= M(m,n); %voisin de droite
      end
  
  else if (n==1)
    P{m+1,n}{1}= M(m,n); %voisin inférieur
    P{m-1,n}{2}= M(m,n); %voisin supérieur
    P{m,n+1}{4}= M(m,n); %voisin de droite
  else if (n==8)
    P{m+1,n}{1}= M(m,n); %voisin inférieur
    P{m-1,n}{2}= M(m,n); %voisin supérieur
    P{m,n-1}{3}= M(m,n); %voisin de gauche
  else
    P{m+1,n}{1}= M(m,n); %voisin inférieur
    P{m-1,n}{2}= M(m,n); %voisin supérieur
    P{m,n-1}{3}= M(m,n); %voisin de gauche
    P{m,n+1}{4}= M(m,n); %voisin de droite
  end
  
end 1









