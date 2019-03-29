clear all

A=[1 -0.2; -0.2 1]; % matrice des coefficients

a=1;


k_bT=0.2; % Temperature arbitraire
Hydrophile=2;
nb_Hydrophile=28;
Hydrophobe=1;
nb_Hydrophobe=64-nb_Hydrophile;


% 
% for i=1:8
%     M(1,i)=2;
%     M(i,1)=2;
%     M(8,i)=2;
%     M(i,8)=2;
% end
% 
% %matrices des voisins cenvention [N,S,E,O]
 P=cell(8);
P{1,1}=[0,1,1,0];
P{1,8}=[0,1,0,1];
P{8,1}=[1,0,1,0];
P{8,8}=[1,0,0,1];
for i=2:7
    P{1,i}=[0,1,1,1];
    P{i,1}=[1,1,1,0];
    P{8,i}=[1,0,1,1];
    P{i,8}=[1,1,0,1];
end
P{2,2}=[1,1,1,1];
P{2,7}=[1,1,1,1];
P{7,2}=[1,1,1,1];
P{7,7}=[1,1,1,1];
for i=3:6
    P{2,i}=[1,1,1,1];
    P{i,2}=[1,1,1,1];
    P{7,i}=[1,1,1,1];
    P{i,7}=[1,1,1,1];
end
for i=3:6
    for j=3:6
        P{i,j}=[1,1,1,1];
    end
end



%matrice des acides aminés
M=ones(9,9);
k=1;
while k<=nb_Hydrophile
    m=randi(8);
    n=randi(8);
    if M(m,n)==1
        M(m,n)=2;
        if (m==1)
            if(n==1)
                 P{m+1,n}(1)= M(m,n); %voisin inférieur
                 P{m,n+1}(4)= M(m,n); %voisin de droite
            elseif(n==8)
                 P{m+1,n}(1)= M(m,n); %voisin inférieur
                 P{m,n-1}(3)= M(m,n); %voisin de gauche
            else
                 P{m+1,n}(1)= M(m,n); %voisin inférieur
                 P{m,n-1}(3)= M(m,n); %voisin de gauche
                 P{m,n+1}(4)= M(m,n); %voisin de droite
            end
         elseif (m==8)
            if (n==1)
                P{m-1,n}(2)= M(m,n); %voisin supérieur
                P{m,n+1}(4)= M(m,n); %voisin de droite
            elseif (n==8)
                P{m-1,n}(2)= M(m,n); %voisin supérieur
                P{m,n-1}(3)= M(m,n); %voisin de gauche
            else
                P{m-1,n}(2)= M(m,n); %voisin supérieur
                P{m,n-1}(3)= M(m,n); %voisin de gauche
                P{m,n+1}(4)= M(m,n); %voisin de droite
            end
  
        elseif (n==1)
            P{m+1,n}(1)= M(m,n); %voisin inférieur
            P{m-1,n}(2)= M(m,n); %voisin supérieur
            P{m,n+1}(4)= M(m,n); %voisin de droite
        elseif (n==8)
            P{m+1,n}(1)= M(m,n); %voisin inférieur
            P{m-1,n}(2)= M(m,n); %voisin supérieur
            P{m,n-1}(3)= M(m,n); %voisin de gauche
            else
            P{m+1,n}(1)= M(m,n); %voisin inférieur
            P{m-1,n}(2)= M(m,n); %voisin supérieur
            P{m,n-1}(3)= M(m,n); %voisin de gauche
            P{m,n+1}(4)= M(m,n); %voisin de droite
        end
        k=k+1;
    end
   
end

En=[];
En(1)=EnergieTotale(P,M,A,a);
iter=1;
eps=1e-6;
arret=0;
Emin=En(1);
Mconfig=zeros(9,9);


 

while (arret==0)
    %tirage aléatoire
    n=randi(8);
    m=randi(8);
    o=randi(8);
    p=randi(8);
    if(M(m,n)~=M(o,p))
        dE=variationE_test(P,M,A,a,m,n,o,p);
        if (dE<0||((dE>0)&&(rand(1)<exp(-(dE)/k_bT))))
            %changement AA
            memoire=0;
            memoire=M(m,n);
            M(m,n)=M(o,p);
            M(o,p)=memoire;
            %changements pour les voisins
            if (m==1)
                if(n==1)
                     P{m+1,n}(1)= M(m,n); %voisin inférieur
                     P{m,n+1}(4)= M(m,n); %voisin de droite
                elseif(n==8)
                     P{m+1,n}(1)= M(m,n); %voisin inférieur
                     P{m,n-1}(3)= M(m,n); %voisin de gauche
                else
                     P{m+1,n}(1)= M(m,n); %voisin inférieur
                     P{m,n-1}(3)= M(m,n); %voisin de gauche
                     P{m,n+1}(4)= M(m,n); %voisin de droite
                end
             elseif (m==8)
                if (n==1)
                    P{m-1,n}(2)= M(m,n); %voisin supérieur
                    P{m,n+1}(4)= M(m,n); %voisin de droite
                elseif (n==8)
                    P{m-1,n}(2)= M(m,n); %voisin supérieur
                    P{m,n-1}(3)= M(m,n); %voisin de gauche
                else
                    P{m-1,n}(2)= M(m,n); %voisin supérieur
                    P{m,n-1}(3)= M(m,n); %voisin de gauche
                    P{m,n+1}(4)= M(m,n); %voisin de droite
                end

            elseif (n==1)
                P{m+1,n}(1)= M(m,n); %voisin inférieur
                P{m-1,n}(2)= M(m,n); %voisin supérieur
                P{m,n+1}(4)= M(m,n); %voisin de droite
            elseif (n==8)
                P{m+1,n}(1)= M(m,n); %voisin inférieur
                P{m-1,n}(2)= M(m,n); %voisin supérieur
                P{m,n-1}(3)= M(m,n); %voisin de gauche
                else
                P{m+1,n}(1)= M(m,n); %voisin inférieur
                P{m-1,n}(2)= M(m,n); %voisin supérieur
                P{m,n-1}(3)= M(m,n); %voisin de gauche
                P{m,n+1}(4)= M(m,n); %voisin de droite
            end


        En(iter+1)=En(iter)+dE;


        else
        En(iter+1)=En(iter);

        end

        iter=iter+1;
        if (iter>201&&(En(iter)<=En(1)))
            if(abs(mean(En(iter-10:iter))-mean(En(iter-200:iter-100)))<eps)
            arret=1;
            end
        end

        if ((En(iter)<Emin))  
            Emin=En(iter);
            Mconfig=M;
        end
    end

%     subplot(1,2,1);pcolor(M); axis square;
% 
%     subplot(1,2,1)
%     refreshdata
%     drawnow
%     pcolor(M); axis square;
%     
%     subplot(1,2,2)
%     refreshdata
%     drawnow
%     plot((1:iter),En); axis square;
%     xlabel('Iteration');
%     ylabel('Potentiel')

end  


subplot(1,2,1);pcolor(M); axis square;

    subplot(1,2,1)
    pcolor(M); axis square;
    subplot(1,2,2)
    plot((1:iter),En); axis square;
    xlabel('Iteration');
    ylabel('Potentiel')


%pcolor(Mconfig);axis square;

