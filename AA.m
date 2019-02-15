
classdef AA
    
    properties
        type % hydrophile ou hydrophobe
        v0 % vitesse initiale
        
        positionx
        positiony
        
    end
    
    methods
        function obj = AA(type,x,y) % initialise les valeurs initiales
            global V;
            obj.type = type;
            obj.v0 = random_maxboltz(V(2),1,V(3)); % V(2) = T & V(3) = m
            obj.positionx = x;
            obj.positiony = y;
        end
        function Ft = force(obj, protein)
            global V epsi;
            
            Ft = zeros(2,1);
            
            for i = 1:length(protein)
                if ~(isequal(protein(i),obj))
                    dx = obj.positionx - protein(i).positionx;
                    dy = obj.positiony - protein(i).positiony;
                    
                    r = sqrt(dx^2+dy^2);
                    
                    theta = atan(dy/dx);
                    
                    % V(9) = r0
                    if  r <= V(9)/sqrt(2)
                        F = LenardJones(V(9),r,epsi(obj.type,protein(i).type),1);
                    else
                        F = 0;
                    end
                    
                    Ft = Ft + [F*cos(theta);F*sin(theta)];
                end
            end
            
        end
        function update(obj,r)
            obj.positionx = obj.positionx + r(1);
            obj.positiony = obj.positionx + r(2);
        end
        function [x,y] = position(obj)
            x = obj.positionx;
            y = obj.positiony;
        end
        
    end
    
end

