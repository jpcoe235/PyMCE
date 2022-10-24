function E = padm( A, p )

if nargin == 1, p = 6; end;
[n,n] = size(A);

% Pade coefficients (1-based instead of 0-based as in the literature)

c(1) = 1;
for k = 1:p
  c(k+1) = c(k)*((p+1-k)/(k*(2*p+1-k)));
end;

% Scaling

s = norm(A,'inf');
if s > 0.5,
  s = max(0,fix(log(s)/log(2))+2);
  A = 2^(-s)*A;
end;

% Horner evaluation of the irreducible fraction (see ref. above)

I = eye(n);
A2 = A*A;
Q = c(p+1)*I;
P = c(p)*I;
odd = 1;
for k = p-1:-1:1,
  if odd == 1,
    Q = Q*A2 + c(k)*I;
  else
    P = P*A2 + c(k)*I;
  end;
  odd = 1-odd;
end;
if odd == 1
  Q = Q*A;
  Q = Q - P;
  E = -(I + 2*(Q\P));
else
  P = P*A;
  Q = Q - P;
  E = I + 2*(Q\P);
end;

% Squaring

for k = 1:s,
  E = E*E;
end;