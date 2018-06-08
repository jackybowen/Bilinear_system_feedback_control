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


function x = mono_upto_next_grlex ( m, n, x )

%*****************************************************************************80
%
%% MONO_UPTO_NEXT_GRLEX: grlex next monomial with total degree up to N.
%
%  Discussion:
%
%    We consider all monomials in an M dimensional space, with total
%    degree up to N.
%
%    For example:
%
%    M = 3
%    N = 3
%
%    #  X(1)  X(2)  X(3)  Degree
%      +------------------------
%    1 |  0     0     0        0
%      |
%    2 |  0     0     1        1
%    3 |  0     1     0        1
%    4 |  1     0     0        1
%      |
%    5 |  0     0     2        2
%    6 |  0     1     1        2
%    7 |  0     2     0        2
%    8 |  1     0     1        2
%    9 |  1     1     0        2
%   10 |  2     0     0        2
%      |
%   11 |  0     0     3        3
%   12 |  0     1     2        3
%   13 |  0     2     1        3
%   14 |  0     3     0        3
%   15 |  1     0     2        3
%   16 |  1     1     1        3
%   17 |  1     2     0        3
%   18 |  2     0     1        3
%   19 |  2     1     0        3
%   20 |  3     0     0        3
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 September 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the maximum degree.
%    0 <= N.
%
%    Input, integer X(M), the current monomial.
%    To start the sequence, set X = [ 0, 0, ..., 0, 0 ].
%
%    Output, integer X(M), the next monomial.
%    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
%
  if ( n < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MONO_UPTO_NEXT_GRLEX - Fatal error!\n' );
    fprintf ( 1, '  N < 0.\n' );
    error ( 'MONO_UPTO_NEXT_GRLEX - Fatal error!' );
  end

  if ( sum ( x(1:m) ) < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MONO_UPTO_NEXT_GRLEX - Fatal error!\n' );
    fprintf ( 1, '  Input X sums to less than 0.\n' );
    error ( 'MONO_UPTO_NEXT_GRLEX - Fatal error!' );
  end

  if ( n < sum ( x(1:m) ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MONO_UPTO_NEXT_GRLEX - Fatal error!\n' );
    fprintf ( 1, '  Input X sums to more than N.\n' );
    error ( 'MONO_UPTO_NEXT_GRLEX - Fatal error!' );
  end

  if ( n == 0 )
    return
  end

  if ( x(1) == n )
    x(1) = 0;
  else
    x = mono_next_grlex ( m, x );
  end

  return
end

function x = mono_next_grlex ( m, x )

%*****************************************************************************80
%
%% MONO_NEXT_GRLEX: grlex next monomial.
%
%  Discussion:
%
%    Example:
%
%    M = 3
%
%    #  X(1)  X(2)  X(3)  Degree
%      +------------------------
%    1 |  0     0     0        0
%      |
%    2 |  0     0     1        1
%    3 |  0     1     0        1
%    4 |  1     0     0        1
%      |
%    5 |  0     0     2        2
%    6 |  0     1     1        2
%    7 |  0     2     0        2
%    8 |  1     0     1        2
%    9 |  1     1     0        2
%   10 |  2     0     0        2
%      |
%   11 |  0     0     3        3
%   12 |  0     1     2        3
%   13 |  0     2     1        3
%   14 |  0     3     0        3
%   15 |  1     0     2        3
%   16 |  1     1     1        3
%   17 |  1     2     0        3
%   18 |  2     0     1        3
%   19 |  2     1     0        3
%   20 |  3     0     0        3
%
%    Thanks to Stefan Klus for pointing out a discrepancy in a previous
%    version of this code, 05 February 2015.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 February 2015
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer X(M), the current monomial.
%    The first item is X = [ 0, 0, ..., 0, 0 ].
%
%    Output, integer X(M), the next monomial.
%

%
%  Ensure that 1 <= M.
%
  if ( m < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MONO_NEXT_GRLEX - Fatal error!' );
    fprintf ( 1, '  M < 1\n' );
    error ( 'MONO_NEXT_GRLEX - Fatal error!' );
  end
%
%  Ensure that 0 <= XC(I).
%
  for i = 1 : m
    if ( x(i) < 0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'MONO_NEXT_GRLEX - Fatal error!' );
      fprintf ( 1, '  X(I) < 0\n' );
      error ( 'MONO_NEXT_GRLEX - Fatal error!' );
    end
  end
%
%  Find I, the index of the rightmost nonzero entry of X.
%
  i = 0;
  for j = m : -1 : 1
    if ( 0 < x(j) )
      i = j;
      break
    end
  end    
%
%  set T = X(I)
%  set X(I) to zero,
%  increase X(I-1) by 1,
%  increment X(M) by T-1.
%
  if ( i == 0 )
    x(m) = 1;
    return
  elseif ( i == 1 )
    t = x(1) + 1;
    im1 = m;
  elseif ( 1 < i )
    t = x(i);
    im1 = i - 1;
  end

  x(i) = 0;
  x(im1) = x(im1) + 1;
  x(m) = x(m) + t - 1;

  return
end
