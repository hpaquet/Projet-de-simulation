function dE=variationE_test(M_voisins,M_AA,M_coeff,a,m,n,o,p)
dE=0;
for k=1:4
    if (M_voisins{m,n}(k)~=0)
        dE=dE-LenardJones(a,a,M_coeff(M_AA(m,n),M_voisins{m,n}(k)),0);
        dE=dE+LenardJones(a,a,M_coeff(M_AA(o,p),M_voisins{m,n}(k)),0);
    end
        
    if (M_voisins{o,p}(k)~=0)
        dE=dE-LenardJones(a,a,M_coeff(M_AA(o,p),M_voisins{o,p}(k)),0);
        dE=dE+LenardJones(a,a,M_coeff(M_AA(m,n),M_voisins{o,p}(k)),0);
        
    end
    
end   
end