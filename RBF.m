function [ Psi ] = RBF( X ,centers,sigma)
%RBF radial basis function
%   Gaussian function used for radial basis function
% X is a scalar or column vector
% output Psi as a row vector of [psi1 psi2 psi3 ...]
% sigma = 0.1;
for i = 1:size(centers,2)
    r = norm(X-centers(:,i));
    Psi(i) = exp(-r^2/sigma^2);% Gaussian kernel
%     Psi(i) = r^2*log(r);% Thin plate spline

end

end
