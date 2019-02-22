function [ Fx, Fy ] = sum_force( protein,  i, V)

    for n = 1:length(protein)
        if n ==i
            continue
            
        elseif n==i-1 || n== i+1
            P = LenardJones(V(6),r,e,d);
            
        else
            P = LenardJones(V(6),r,e,d);
            
        end
        
        

    end
    


end

