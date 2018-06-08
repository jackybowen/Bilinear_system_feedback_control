function [ C_vec ] = coeffs_nD_PDE(P, Psi, n)
%COEFFS_nD 
%~~~~~~~~~~Extract coefficients of Polynomial term corresponding to 
% given nD monomial basis.
%   Input: P - - - - - - The polynomial term to be expanded
%          Monom - - - - The monomials in n variables by natural ordering
%   Output:
%          Cxy - - - - - The coefficient row vector of polynomial in the 
%   given monomial basis space
N = length(Psi);

x = sym('y',[n,1]);
C_vec = zeros(1,N);
[Cxy,Txy] = coeffs(P,x);
if ~isempty(Cxy)
    Cxy = reshape(Cxy,1,[]);
    Txy = reshape(Txy,1,[]);
    for i = 1:length(Txy)
        C_vec(Psi==Txy(i)) = Cxy(i);
    end
end

end

