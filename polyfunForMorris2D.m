function Y = polyfunForMorris2D(X)
% evaluate polynomial
% x is N X D

% simply quadratic
Y = 3*X(:,1) + X(:,2).^2;

