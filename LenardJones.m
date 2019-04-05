
function P = LenardJones(a,r,e,d)

r0 = (3-sqrt(2))^(-1/6)*a;

switch d
    case 0 % Potentiel
        
        P = ( (r0./r).^12 - e.*(r0./r).^6 ).*(1-3/(2*a^2).*r.^2 + 1/(sqrt(2)*a^3).*r.^3);
        
    case 1 % force
        
        %r = double(r);
        
        rm7 = 1./(r.^7);
        rm13 = 1./(r.^13);
        
        P = -( ( e*6*r0^6.*rm7 - 12*r0^12.*rm13).*(1-0.5*3*a^(-2)*r.^2 + 1/(2^(1/2)*a^3)*r.^3) + ...
            ( r0.^12.*rm13.*r - e*(r0).^6.*r.*rm7 ).*(-3*a^(-2)*r + 3/(2^(1/2)*a^3)*r.^2));
        
    case 2 % constante de ressord
        
        r = a;
        
        P = (( (156*r0^12)./(r.^14) - e*(42*r0^6)./(r.^8) ).*(1-3/(2*a^2)*r.^2 + 1/(sqrt(2)*a^3)*r.^3) + ...
            2*( e*((6*r0^6)./(r.^7)) - (12*r0^12)./(r.^13)).*(-3/(a^2)*r + 3/(sqrt(2)*a^3)*r.^2) + ...
            (-3/a^2 + 6/(sqrt(2)*a^3)*r).*( (r0./r).^12 - e*(r0./r).^6 ));
        
end

end
