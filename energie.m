% Calcul l'énergie total du système
function [K E] = energie(x,y,a,N,t,m,dt)
    
    global k VOI

    di = eye(N)>0;

    dx = (x(t,:).*ones(N,N))-x(t,:)';dy = (y(t,:).*ones(N,N))-y(t,:)'; % Distance x et y entre les AAs
    r = sqrt(dx.^2+dy.^2); % Distance entre les AAs
    delta = r-a; % Distance du point d'équilibre des AAs
    
    U = VOI.*0.5.*k.*delta.^2;
    U(di) = 0;
    U = 0.5*sum(sum(U,2));
    
    vx = (3*x(t,:)-4*x(t-1,:)+x(t-2,:))/(2*dt);
    vy = (3*y(t,:)-4*y(t-1,:)+y(t-2,:))/(2*dt);
    v2 = vx.^2+vy.^2;
    
    K = 0.5*m*sum(v2);
    
    E = K + U;

end