function [Psi] = RBF(X, centers)
%RBF radial basis function
%   Gaussian function used for radial basis function
%Input:
%       X -------------------- input data(n x 1)
%    centers ----------------- N centers data(n x N)
%Output:
%       Psi is a row vector of [psi_1 psi_2 psi_3 ... psi_N]
sigma = 0.1;
for i = 1:size(centers,2)
    r = norm(X-centers(:,i));
    Psi(i) = exp(-r^2/sigma^2); % Gaussian kernel
%     Psi(i) = r^2*log(r);% Thin plate spline

end

end
