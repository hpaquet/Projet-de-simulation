function [ Fx, Fy ] = sum_force( protein, type, i, V,t,e)

    l = size(protein,2);
    
    Fx = 0;
    Fy = 0;

    for n = 1:l/2
        
        if n ==i
            continue
            
        elseif n==i-1 || n== i+1
            
            dx = (protein(t,2*i-1) - protein(t,2*n-1));
            dy = protein(t,2*n) - protein(t,2*i);
            
            r = sqrt(dx^2+dy^2);
            
            theta = atan(dy/dx);
            
            k = 1;%LenardJones(V(6),r,e(type(i),type(n)),2)*10;
            
            F = -k*(dx-V(6));
            
            Fx = Fx + F*cos(theta);
            Fy = Fy + F*sin(theta);
            
%         else
%             dx = protein(t,2*n-1) - protein(t,2*i-1);
%             dy = protein(t,2*n) - protein(t,2*i);
%             
%             r = sqrt(dx^2+dy^2);
%             
%             theta = atan(dy/dx);
%             
%             F = -LenardJones(V(6),r,e(type(i),type(n)),1);
%             
%             Fx = Fx + F*cos(theta);
%             Fy = Fy + F*sin(theta);
            
        end
        
        

    end
    


end

