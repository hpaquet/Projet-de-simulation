function P=changement_voisins(P,M,m,n)


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

end