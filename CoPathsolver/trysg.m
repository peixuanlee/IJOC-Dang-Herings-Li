clear sum;
global s n m  m0 n1 m1 pm0 pm1 MU delta0 v0 timep nn fes


%Inputs: s, n, m, u, tpm

delta0 = 0.95;
%Discount factor


n1 = s*n;
m0 = ones(1, s)*m*ones(n, 1);
m1 = m*ones(n, 1);

MU = zeros(m0, n1);
%Matrix corresponding to mu
q0 = 0;
q1 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        MU(h0+1:h0+m(q, i), q1+i) = 1;
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
    q1 = q1 + n;
end


pm0 = zeros(s, 1);
for q = 1:s
    pm0(q) = prod(m(q, :));
end
pm1 = zeros(s, 1);
for h = 1:s
    pm1(h) = pm0(h)/min(m(h,:));
end


v0 = zeros(m0, 1);
%Starting point
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        wa = rand(m(q, i), 1);       
        v0(h0+1:h0+m(q, i)) = wa/sum(wa);    
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end




tic;
xx = pathmcp([v0;rand(n1,1)],[zeros(m0,1);-inf*ones(n1,1)],inf*ones(m0+n1,1),'tryexm');
%xx = pathmcp(rand(m0+n1,1),[zeros(m0,1);-inf*ones(n1,1)],inf*ones(m0+n1,1),'tryexm')


fprintf('The Dimension of Variables                 = %15.2f\n\n', length(xx))
f = tryexm(xx,0);

resid = f'*xx;
timep = toc;

sum = sum(xx(1:m(1, 1)));
if sum >= 1.001|| sum <= 0.999 || fes == 0
    nn = 0; timep = 0;
    display 'the path fails to find a solution';
else   
    if resid >1.0e-6 
        display 'the path finds a bad solution';
    end
    nn = 1;
    fprintf('Computational Time of PathSolver                 = %15.2f\n\n', timep);
end
fprintf('Stationary Equilibrium');
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        fprintf('\n');
        fprintf('x =\n');
        fprintf('%15.12f %15.12f %15.12f %15.12f %15.12f', xx(h0+1:h0+m(q, i)));
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end
fprintf('\n\n');

homo;

fprintf('\n');

fprintf('___________________________________________________________________________\n');


