
function [x,y] = simulation(dt,tmax,T0,N,m,a,THERMO,v0,x,y,fig_on)
%% Paramètre graphique
if fig_on
    
% iteration par image
IPA = 100;
    
% Dimension du graphique
xx = N*1e-9;
yy = N*1e-9;

figure('Name','Dynamique moléculaire','Position', [ 50 50 1200 600 ],'NumberTitle','off');
f1 = subplot(2,3,[1 2 4 5]);
h1=plot(0,0,'MarkerSize',50,'Marker','.','LineWidth',3);
axis off
axis(f1,[-xx xx -yy yy])
set(gca,'nextplot','replacechildren')
subplot(2,3,3)
h2=animatedline;
title('Température en fonction du temps')
ylabel('Température (K)')
xlabel('Temps (as)')
subplot(2,3,6)
h3=animatedline;
title('Énergie total en fonction du temps')
ylabel('Énergie (J)')
xlabel('Temps (as)')

end

%% Première itération (t = 1)

C3 = dt^2*force(x(1,:),y(1,:),a,N)/(2*m);

x(2,:) = x(1,:) + C3(1,:) + dt*v0(1,:);
y(2,:) = y(1,:) + C3(2,:) + dt*v0(2,:);


%% Boucle principale (pour t>1)

t = 2;

while(t < tmax)
    
    % Mise à jour des positions
    
    C3 = dt^2*force(x(t,:),y(t,:),a,N)/m;
    
    x(t+1,:) = 2*x(t,:) + C3(1,:) - x(t-1,:);
    y(t+1,:) = 2*y(t,:) + C3(2,:) - y(t-1,:);
    
    t = t+1;
    
    % Température du système
    T = temperature(x,y,N,dt,m,t);
   
    % Énergie du système
    E = energie(x,y,a,N,t,m,dt);
    
    %  Mise à jour des graphique à tout les 10 images
    if fig_on && mod(t,IPA) == 0 
        
        [xs,ys] = mc(x,y,N,t); % calcul du centre de masse
        
        addpoints(h2,t*dt/1e-18,T); % maj du graphe de la temperature
        
        addpoints(h3,t*dt/1e-18,E); % maj du graphe de la temperature
        
        set(h1,'XData',x(t,:),'YData',y(t,:)); % maj des positions
        axis(f1,[-xx+xs xx+xs -yy+ys yy+ys]); % maj des axes
        
        drawnow;
    end
    
    %  Thermostat
    if THERMO && (T > 10*T0 || T< T0/10 ) && t>3
        
        alp2 = sqrt(T0/T); % coefficient alpha
        
        % Valeur des vitesses
        vx = (3*x(t,:) - 4*x(t-2,:)  + x(t-3,:))./(2*dt);
        vy = (3*y(t,:) - 4*y(t-2,:)  + y(t-3,:))./(2*dt);
        
        % Ajustement des vitesses selon alpha
        v0 = alp2 * [vx; vy];
        
        % Première itération de la simulation pour repartir la simulation
        C3 = dt^2*force(x(1,:),y(1,:),a,N)/(2*m);
        
        x(t+1,:) = x(t,:) + C3(1,:) + dt*v0(1,:);
        y(t+1,:) = y(t,:) + C3(2,:) + dt*v0(2,:);
        
        t = t+1;
        
    end
    
    % initialise 1000 nouveau espace mémoire dans les variable de position
    if mod(t,1000) == 0
        x = [x;zeros(1000,N)];
        y = [y;zeros(1000,N)];
    end
    
end

close;

end




