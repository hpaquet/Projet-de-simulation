function protein = MD_simulation(protein, type, duration, v0)

global V epsi;

N = size(protein,2)/2;

h=plot(0,0,'MarkerSize',100,'Marker','.');

set(gca,'nextplot','replacechildren');


%% première itération

for n = 1:N % ensemble des AA
    [Fx, Fy] = sum_force( protein, type, n, V, 1, epsi);
    
    protein(2,2*n-1) = 2*V(9)/V(4)*Fx + 2*V(9)*sqrt(2*V(5))*V(7) + protein(1,2*n-1) - V(9)*v0(n);
    protein(2,2*n) = 2*V(9)/V(4)*Fy + 2*V(9)*sqrt(2*V(5))*V(7) + protein(1,2*n) - V(9)*v0(n);
end

%% boucle principal

for t = 3:duration % variation dans le temps
    
    for n = 1:N % ensemble des AA
        
        [Fx, Fy] = sum_force( protein, type, n, V, t-1, epsi);
        
        protein(t,2*n-1) = 2*V(9)/V(4)*Fx + 2*V(9)*sqrt(2*V(5))*V(7) + protein(t-1,2*n-1);
        protein(t,2*n) = 2*V(9)/V(4)*Fy + 2*V(9)*sqrt(2*V(5))*V(7) + protein(t-1,2*n);
            
        
    end
    
    set(h,'XData',[protein(t,1) protein(t,3)],'YData',[protein(t,2) protein(t,4)]);
    drawnow;
    pause(0.1)
    
end

end

