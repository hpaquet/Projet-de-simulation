clc;close all;clear all;
%%
% Fonctionne parfaitement
%
%
%
%%

m = 1;
k = 10;
d = 0.005;
a = 4;
g = 0;
xi = 0;
D = 0;
tf = 100;
j = 1;

x_sol = [];

%%

for i = [1 1/2 1/4 1/8 1/16 1/32 1/64 1/128 1/256 1/512]
    
    dt = i*d;
    N = tf/i;
    
    for n = 1:2
        x(1,n) = n;
        v0(1,n) = (-1)^(n-1);
    end
    
    for n = 1:2
    C3 = dt^2*force(x,n,k,a);

    x(2,:) = x(1,:) + C3(1,:) + dt*v0(1,:);
    end
    
    
    for t = 2:N
        
        for n =1:2
        C3 = dt^2*force(x,n,k,a)/m;

        x(t+1,:) = 2*x(t,:) + C3(1,:) - x(t-1,:);
            
        end
        
    end
    
    x_sol(end+1,:) = sum(x,1)/size(x,1);
    
end

err = [];

for k = 1:length(x_sol)-1
    err(k,:) = abs(x_sol(k,:) - x_sol(k+1,:));
end

N = [1 1/2 1/4 1/8 1/16 1/32 1/64 1/128 1/256 1/512]';
 
loglog(N(2:end),err)
polyfit(log(N(2:end)),log(err(:,1)),1)
polyfit(log(N(2:end)),log(err(:,2)),1)


%%

function F = force(x,n,k,a)

F = 0;

if n-1 >= 1
    F = F - k*(x(n)-x(n-1)-a);
end
if n+1 <= length(x)
    F = F - k*(x(n)-x(n+1)+a);
end

end
