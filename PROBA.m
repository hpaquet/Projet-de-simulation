clear all
A=[1 -0.5; -0.5 -1]; % matrice des coefficients

a=1;


k_bT=0.035; % Temperature arbitraire

nb_Hydrophile=32;

nb_Hydrophobe=64-nb_Hydrophile;


% 
% for i=1:8
%     M(1,i)=2;
%     M(i,1)=2;
%     M(8,i)=2;
%     M(i,8)=2;
% end
% 
% %matrices des voisins cenvention [N,S,E,O]
 P=cell(8);
P{1,1}=[0,1,1,0];
P{1,8}=[0,1,0,1];
P{8,1}=[1,0,1,0];
P{8,8}=[1,0,0,1];
for i=2:7
    P{1,i}=[0,1,1,1];
    P{i,1}=[1,1,1,0];
    P{8,i}=[1,0,1,1];
    P{i,8}=[1,1,0,1];
end
P{2,2}=[1,1,1,1];
P{2,7}=[1,1,1,1];
P{7,2}=[1,1,1,1];
P{7,7}=[1,1,1,1];
for i=3:6
    P{2,i}=[1,1,1,1];
    P{i,2}=[1,1,1,1];
    P{7,i}=[1,1,1,1];
    P{i,7}=[1,1,1,1];
end
for i=3:6
    for j=3:6
        P{i,j}=[1,1,1,1];
    end
end

M={};
for i=1:1000
    M{i}=simu(k_bT,nb_Hydrophile,nb_Hydrophobe,P,A);
end



%%
RES=ones(9,9);
for i=1:8
    for j=1:8
        S=0;
        for k=1:1000
            S=S+M{k}(i,j)-1; 
        end
        RES(i,j)=S/1000;
    end
end
 
pcolor(RES)
            
            
        