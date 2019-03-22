function E=EnergieTotale(M_voisins,M_AA,M_coeff,a)
E=0;
for i=1:8
    for j=1:8
        for k=1:4
            if (M_voisins{i,j}(k)~=0)
            E=E+LenardJones(a,a,M_coeff(M_AA(i,j),M_voisins{i,j}(k)),0);
            end
        end
    end
end
end