function Mconfig=simu(k_bT,nb_Hydrophile,nb_Hydrophobe,P,A)

rng('shuffle')


a=1;
Hydrophile=2;
Hydrophobe=1;

%matrice des acides aminés
M=ones(9,9);
k=1;
while k<=nb_Hydrophile
    m=randi(8);
    n=randi(8);
    if M(m,n)==1
        M(m,n)=2;
        P=changement_voisins(P,M,m,n);
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
            P=changement_voisins(P,M,m,n);
            P=changement_voisins(P,M,o,p);

        En(iter+1)=En(iter)+dE;


        else
        En(iter+1)=En(iter);

        end

        iter=iter+1;
        if (iter>201&&(En(iter)<=En(1)))
        %if  (iter>=100)
            if(abs(mean(En(iter-100:iter))-mean(En(iter-200:iter-100)))<eps)
            arret=1;
            end
        end

        if ((En(iter)<Emin))  
            Emin=En(iter);
            Mconfig=M;
        end
    end





end