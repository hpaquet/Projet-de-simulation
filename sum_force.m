function [ Fx, Fy ] = sum_force( protein, type, i, V,t,e)

    l = length(protein);
    
    Fx = 0;
    Fy = 0;

    for n = 1:l
        
        if n ==i
            continue
            
        elseif n==i-1 || n== i+1
            
            dx = protein(2*n-1,t) - protein(2*i-1,t);
            dy = protein(2*n,t) - protein(2*i,t);
            
            r = sqrt(dx^2+dy^2);
            
            theta = atan(dy/dx);
            
            k = LenardJones(V(6),r,e(type(i),type(n)),2)*10;
            
            F = -0.5*k*r^2;
            
            Fx = Fx + F*cos(theta);
            Fy = Fy + F*sin(theta);
            
        else
            dx = protein(2*n-1,t) - protein(2*i-1,t);
            dy = protein(2*n,t) - protein(2*i,t);
            
            r = sqrt(dx^2+dy^2);
            
            theta = atan(dy/dx);
            
            F = -LenardJones(V(6),r,e(type(i),type(n)),1);
            
            Fx = Fx + F*cos(theta);
            Fy = Fy + F*sin(theta);
            
        end
        
        

    end
    


end

