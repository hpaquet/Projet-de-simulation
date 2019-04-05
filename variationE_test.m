function dE=variationE_test(P,M,M_coeff,a,m,n,o,p)
dE=0;
for k=1:4
    if (P{m,n}(k)~=0)
        dE=dE-LenardJones(a,a,M_coeff(M(m,n),P{m,n}(k)),0);
        %u=LenardJones(a,a,M_coeff(M(m,n),P{m,n}(k)),0);
    end
        
    if (P{o,p}(k)~=0)
        dE=dE-LenardJones(a,a,M_coeff(M(o,p),P{o,p}(k)),0);
        %LenardJones(a,a,M_coeff(M(o,p),P{o,p}(k)),0);
        
    end
end

memoire=0;
memoire=M(m,n);
M(m,n)=M(o,p);
M(o,p)=memoire;


P=changement_voisins(P,M,m,n);
P=changement_voisins(P,M,o,p);
M(o,p)=M(m,n);
M(m,n)=memoire;

for k=1:4  
    if (P{m,n}(k)~=0)
        dE=dE+LenardJones(a,a,M_coeff(M(o,p),P{m,n}(k)),0);
        %LenardJones(a,a,M_coeff(M(o,p),P{m,n}(k)),0);
    end
     if (P{o,p}(k)~=0)
        dE=dE+LenardJones(a,a,M_coeff(M(m,n),P{o,p}(k)),0);
        %LenardJones(a,a,M_coeff(M(m,n),P{o,p}(k)),0);
        
     end
end   
P=changement_voisins(P,M,m,n);
P=changement_voisins(P,M,o,p);

end