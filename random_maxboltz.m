function [ v ] = random_maxboltz(T, N, m)
%
% Fonction qui génère une distribution de vitesse de N particule de masse m 
% selon la distribution de Maxwell-Boltzmann à une température T
%

v = [];
k = 1.38e-23;
a = sqrt(k*T/m); % vitesse la plus probable

for i = 1:N
    p = rand(); % nombre aléatoire entre 0 et 1
    
    %Fonction de répartition de la distribution de Maxwell-Boltzmann
    func = @(x) erf(x/(sqrt(2)*a))-sqrt(2/pi)*(x.*exp(-(x.^2)/(2*a^2)))/(a)-p;
    
    v(end+1) = fzero(func,a); % vitesse associé à la probabilité
end

end
