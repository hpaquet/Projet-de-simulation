function dE=variationE(M_voisins,M_AA,M_coeff,a,m,n)
dE=0;
for k=1:4
    if (M_voisins{m,n}(k)~=0)
            dE=dE-LenardJones(a,a,M_coeff(M_AA(m,n),M_voisins{m,n}(k)),0);
    
        if (M_AA(m,n)==1)
            dE=dE+LenardJones(a,a,M_coeff(2,M_voisins{m,n}(k)),0);
    
        else
            dE=dE+LenardJones(a,a,M_coeff(1,M_voisins{m,n}(k)),0);
        end
    end
    
end   
end