function [  ] = MD_simulation(protein, duration)

global V;

N = length(protein);


%% premi�re it�ration

for n = 1:N % ensemble des AA
    F = Force();
    protein(2*n-1,2) = 2*V(10)/V(4)*F(1) + 2*V(10)*sqrt(2*V(0.2))*V(7) + protein(2*n-1,1) - V(10)*vo;
    protein(2*n,2) = 2*V(10)/V(4)*F(2) + 2*V(10)*sqrt(2*V(0.2))*V(7) + protein(2*n-1,1) - V(10)*vo;
end

%% boucle principal

for t = 2:duration % variation dans le temps
    
    
    for n = 1:N % ensemble des AA
        
        
    end
end


end

