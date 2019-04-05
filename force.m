% Calcul la force entre les AAs
function F = force(x,y,a,N)

global VOI ALT k

di = eye(N)>0;
 % Constante de rapel

dx = (x.*ones(N,N))-x';dy = (y.*ones(N,N))-y'; % Distance x et y entre les AAs
r = sqrt(dx.^2+dy.^2); % Distance entre les AAs
delta = r-a; % Distance du point d'équilibre des AAs

% Force entre les voisins (liaisons covalente)
Fx = VOI.*k.*delta.*dx./r; Fy = VOI.*k.*delta.*dy./r;
Fx(di) = 0;Fy(di) = 0;
Fx = sum(Fx,2);Fy = sum(Fy,2);
F = [Fx';Fy'];

% Force entre les autres
Fx = ALT.*LenardJones(a,delta,ALT,1).*dx./r;Fy = ALT.*LenardJones(a,delta,ALT,1).*dy./r;
Fx(ALT<1) = 0;Fy(ALT<1) = 0;
Fx = sum(Fx,2);Fy = sum(Fy,2);
F = F+[Fx';Fy'];

end