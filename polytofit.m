function [b, x] = polytolfit(f, eps)
  % given a function f and number eps, the function finds a polynomial
  % p interpolating f, trying to minimize its degree while guaranteeing
  % that |p(x)-f(x)| <= eps for -1 <= x <= 1. It is assumed that f is
  % differentiable and its even-degree derivatives do not exceed 6/5,
  % while odd ones do not exceed 8/5. The result is returned as two
  % vectors, b and x, with b representing p's coefficients in a
  % Newton's polynomial base associated with the numbers contained
  % in x. 
  
  n = find_minimal_n(eps);
  % using a theoretical bound for max |f(x)-p(x)| for interpolation
  % based on Chebyshev nodes. Finding smallest possible n that
  % guarantees error smaller than eps
  
  nodes = get_n_chebyshev_nodes(n);
  % chebyshev nodes are best for a guaranteed minimization of |f(x)-p(x)|
  % with a given number of interpolation nodes
  
  values = arrayfun(f, nodes);
  b = newton_interpolate(nodes, values);
  % having found the nodes, we calculate p's coefficients
  
  x = nodes;
end

function b = get_bound_for_n(n)
  if mod(n, 2) == 1
    b = 1.6 / (factorial(n+1));
  else
    b = 1.2 / (factorial(n+1));
  end
end

function n = find_minimal_n(eps)
  n = 0;
  while get_bound_for_n(n) > eps
    n = n+1;
  end
end

function nodes = get_n_chebyshev_nodes(n)
  for k = 1:n
    nodes(k) = cos(pi * (2*k-1)/(2*n));
  end
  nodes = nodes;
end

function coefficients = newton_interpolate(nodes, values)
  % divided differences method
  nodes_no = size(nodes);
  nodes_no = nodes_no(1, 2);
  table = zeros(nodes_no);
  for i = 1:nodes_no
    table(i, 1) = values(i);
  end
  for col = 2:nodes_no
    for row = col:nodes_no
      table(row, col) = (table(row, col-1)-table(row-1, col-1)) / (nodes(row)-nodes(row-col+1));
    end
  end
  for i = 1:nodes_no
    coefficients(i) = table(i, i);
  end
end


