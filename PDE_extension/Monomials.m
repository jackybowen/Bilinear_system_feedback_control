function [ Psi ] = Monomials( X )
% MONOMIALS
%   Generate evaluation at monomial basis functions of most degree D
%   The monomials are sorted in grlex
global Deg
% D = 6; % Most degree D


n = length(X);
X = reshape(X,1,[]);
expon = zeros(1,n);
Psi = [1];
for i = 1:nchoosek(n+Deg,Deg)-1
    expon = mono_upto_next_grlex(n,Deg,expon);
    Psi = [Psi prod(X.^expon)];
end

end