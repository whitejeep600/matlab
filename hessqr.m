
function [hQ, R] = hessqr(A)
  % The function computes the QR decomposition of a hessian matrix A
  % of size (n+1)x(n). Each element under the diagonal is zeroed out
  % with the use of a Givens rotation. The matrix Q is the product of
  % all the rotations used. It is returned in compressed form of size (n)x(2),
  % with the i-th row containing the c and s coefficients of the i-th rotation.
  % The R matrix is the result of applyingg the rotations to A.
  n = size(A, 2);
  if n+1 ~= size(A, 1);
    error("Invalid size of the input matrix")
  end
  if ~ istriu(A(2:(size(A, 1)), :)) 
    error("input matrix not hessian")
  end
  hQ = zeros(n, 2);
  temp1 = zeros(n, 1);
  temp2 = zeros(n, 1);
  for i=1:n;
    diagonal = A(i, i);
    subdiagonal = A(i+1, i);
    r = sqrt(diagonal^2+subdiagonal^2);
    if r == 0
      % identity rotation - no need to zero anything out
      hQ(i, 1) = 1;
      hQ(i, 2) = 0;
      continue
    end
    c = diagonal / r;
    s = -1 * subdiagonal / r;
    hQ(i, 1) = c;
    hQ(i, 2) = s;
    % applying the rotation to the affected rows, storing the results in
    % temporary vectors
    temp2 = (s .* A(i, :)) + (c .* A(i+1, :));
    temp1 = (c .* A(i, :)) - (s .* A(i+1, :));
    A(i,:)   = temp1; % :)
    A(i+1,:) = temp2;
  end
  R = A;
  R = R - tril(R, -2) % might as well, since whatever's left there is approximation error
end


function Y = hessmhq(hQ, A)
  % multiplies the matrix hQ, in the format as described in the hessqr function,
  % by the matrix A
  n = size(hQ, 1);
  if size(hQ, 2) ~= 2
    error("invalid format of the hQ matrix")
  end
  if size(A, 1) ~= n+1
    error("invalid size of the A matrix")
  end
  m = size(A, 2);
  for k = 1:n;
    i = n + 1 - k;
    
    % applying Givens rotations to the affected rows
    c = hQ(i, 1);
    s = hQ(i, 2);
    temp1 = (c .* A(i, :)) + (s .* A(i+1, :));
    temp2 = ((-1*s) .* A(i, :)) + (c .* A(i+1, :));
    A(i,:)   = temp1;
    A(i+1,:) = temp2;
  end
  Y = A;
end

function Y = hessmhqinv(hQ, A)
  % multiplies the inverse of the matrix hQ, in the format as described in the
  % hessqr function, by the matrix A
  n = size(hQ, 1);
  m = size(A, 2);
  if size(hQ, 2) ~= 2
    error("invalid format of the hQ matrix")
  end
  if size(A, 1) ~= n+1
    error("invalid size of the A matrix")
  end
  
  % givens rotations are orthogonal, so the inverse of hQ is the product of the
  % conjugate transposes of its rotations in reverse order, and the conjugate
  % transpose of a givens rotation matrix with coefficients c, s is the same
  % matrix with coefficients c, -s
  hQ(:, 1) = flipud(hQ(:, 1));
  hQ(:, 2) = (-1) .* flipud(hQ(:, 2));
  
  for k = 1:n;
    i = n + 1 - k;
    c = hQ(i, 1);
    s = hQ(i, 2);
    % applying Givens rotations to the affected rows
    temp1 = (c .* A(k, :)) + (s .* A(k+1, :));
    temp2 = ((-1*s) .* A(k, :)) + (c .* A(k+1, :));
    A(k,:)   = temp1;
    A(k+1,:) = temp2;
  end
  Y = A;
end

