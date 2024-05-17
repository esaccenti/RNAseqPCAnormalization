function [PF, comb] = TuckerPhiMat(L)

%Calculated pairwise congruence of Loadings vecotors in L

% Get number of loadings
dim = size(L,2);

% Find all possible pairs
comb = nchoosek(1:dim,2);

% Loop on all possible pairs and calculate Tucker's Phi
for i = 1 : size(comb,1)

    l1 = comb(i,1);
    l2 = comb(i,2);

    PF(i,1) = tf(L(:,l1),L(:,l2));

end

end

%%
function phi = tf(x,y)

 phi = sum( x.*y )/sqrt( sum(x.^2)*sum(y.^2) );

end
