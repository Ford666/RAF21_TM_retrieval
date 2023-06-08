function [ invA, inv_S2 ] = Tikinv( A, p, varargin )
% A the original matrix to be inverted
% p controls the Tikhonov parameter

[~, S2, V] = svd(A'*A);                                                    
S = sqrt(diag(S2));
maxsin = max(S);

if nargin == 1
    p = 0.05;
end

Tiklambda = p*maxsin;

inv_S2 = 1./( S.^2 + Tiklambda^2 );                                                                     
invA = V * diag(inv_S2) * V' * A';                                          
% A' already contain one Sigma, so inv_S2 should cancel one and leave one inverse

end