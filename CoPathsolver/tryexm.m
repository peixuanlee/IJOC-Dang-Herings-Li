function [F, J, domerr] = tryexm(z, jacflag)

%Example in Herings's paper
global m0 n1



% initialize
%z = z(:);
F = sparse(tryfunc(z));
A = sparse(trycons(z));
b = zeros(m0+n1,1);
J = [];
domerr = 0;

if (jacflag)
  J = sparse(tryjac(z));
end

return

