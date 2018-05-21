function [ Psi ] = Monomials( X )
% MONOMIALS
%   Generate evaluation at monomial basis functions of most degree D
%   The monomials are sorted in grlex
global Deg
% D = 6; % Most degree D
n = length(X);
Monom = repmat(reshape(X,[],1),1,Deg+1); 
Monom(:, 1) = 1;
Monom = cumprod(Monom,2);
Monom_cell = mat2cell(Monom,ones(1,n));
Monom_grid = cell(1,numel(Monom_cell));
[Monom_grid{:}] = ndgrid(Monom_cell{:});
Monom = Monom_grid{1};
if n>1
    for i = 2:1:n
        Monom = Monom.*Monom_grid{i};
    end
end
Psi = cell(1,Deg+1);
sub = cell(1,n);

for idx = 1:numel(Monom)
    [sub{:}] = ind2sub((Deg+1)*ones(1,n),idx);
    sum_sub = sum(cell2mat(sub))-n;
    if  sum_sub <= Deg
        if isempty(Psi{sum_sub+1})
            Psi{sum_sub+1} = Monom(idx);
        else
            Psi{sum_sub+1} = [Psi{sum_sub+1}, Monom(idx)];
        end
    end
end
Psi = [Psi{:}];
end