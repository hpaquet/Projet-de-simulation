function [  ] = MD_simulation(protein, duration, movie, fps)

global V;

N = length(protein);

t_max = duration*fps;

X(1,:) = [protein().positionx];
Y(1,:) = [protein().positiony];

%% première itération

for n = 1:N % ensemble des AA
    r = (2*V(10)/V(4))*protein(n).force(protein) + 2*V(10)*sqrt(2*V(5))*V(7)+protein(n).position() - V(10)*protein(n).v0;
    protein(n).positionx = protein(n).positionx + r(1);
    protein(n).positiony = protein(n).positiony + r(2);
end

%% boucle principal

for t = 2:t_max % variation dans le temps
    
    X(t,:) = [protein().positionx];
    Y(t,:) = [protein().positiony];
    
    for n = 1:N % ensemble des AA
        
        r = (2*V(10)/V(4))*protein(n).force(protein) + 2*V(10)*sqrt(2*V(5))*V(7)+protein(n).position() + [X(t-1,n);Y(t-1,n)];
        protein(n).positionx = protein(n).positionx + r(1);
        protein(n).positiony = protein(n).positiony + r(2);
        
    end
end

%% Création d'un fichier vidéo externe
if movie==true
    video = VideoWriter('MD_simulation.avi');
    video.FrameRate = fps;
    open(video)
    writeVideo(video,F);
end

end

