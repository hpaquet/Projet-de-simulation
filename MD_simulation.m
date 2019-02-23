function protein = MD_simulation(protein, type, duration, v0)

global V epsi;

N = size(protein,2)/2;
r = 3;
h=plot(0,0,'MarkerSize',100,'Marker','.');
axis([-1 3 -1 1]);
set(gca,'nextplot','replacechildren');


%% première itération

% for n = 1:N % ensemble des AA
%     [Fx, Fy] = sum_force( protein, type, n, V, 1, epsi);
%     
%     protein(2,2*n-1) = 2*V(9)/V(4)*Fx + 2*V(9)*sqrt(2*V(5))*V(7) + protein(1,2*n-1) - V(9)*v0(n);
%     protein(2,2*n) = 2*V(9)/V(4)*Fy + 2*V(9)*sqrt(2*V(5))*V(7) + protein(1,2*n) - V(9)*v0(n);
% end
protein(2,[1 3]) = [1 -1]*2*V(9)/V(4)*( - ( protein(1,1)-protein(1,3)-V(6) ) ) + 2*V(9)*sqrt(2*V(5))*V(7) + protein(1,[1 3])- V(9)*v0;


%% boucle principal

% V = [N T m g D a xi rayon dt];

for t = 3:duration % variation dans le temps
    
    %for n = 1:N % ensemble des AA
        
        %[Fx, Fy] = sum_force( protein, type, n, V, t-1, epsi);
        
        %protein(t,2*n-1) = 1/(2*V(3)-V(4)*V(9))*(4*V(3)*protein(t,2*n-1))
        %protein(t,2*n-1) = 2*V(9)/V(4)*Fx + 2*V(9)*sqrt(2*V(5))*V(7) + protein(t-2,2*n-1);
        %protein(t,2*n) = 2*V(9)/V(4)*Fy + 2*V(9)*sqrt(2*V(5))*V(7) + protein(t-2,2*n);
            
        protein(t,[1 3]) = [1 -1]*2*V(9)/V(4)*( - ( protein(t,1)-protein(t,3)-V(6) ) ) + 2*V(9)*sqrt(2*V(5))*V(7) + protein(t-2,[1 3]);
        
        
    %end
    
    set(h,'XData',[protein(t,1) protein(t,3)],'YData',[0,0]);
    drawnow;
    pause(0.1)
%     
end

end

