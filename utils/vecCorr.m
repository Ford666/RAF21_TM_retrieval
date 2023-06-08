function pcc = vecCorr(x, y, dim)
% calculate pcc between the rows or columns of two matrix

[rows, cols] = size(x);
if nargin==2 || (nargin>2 && dim==1)
    xi = x; yi = y;
    xi = xi - sum(xi, 1)./rows; yi = yi - sum(yi, 1)./rows;
    pcc = sum(xi.*yi, 1) ./ (vecnorm(xi, 2, 1) .* vecnorm(yi, 2, 1));

elseif (nargin>2 && dim==2)
    xi = x; yi = y;
    xi = xi - sum(xi, 2)./cols; yi = yi - sum(yi, dim)./rows;
    pcc = sum(xi.*yi, 2) ./ (vecnorm(xi, 2, 2) .* vecnorm(yi, 2, 2));
end