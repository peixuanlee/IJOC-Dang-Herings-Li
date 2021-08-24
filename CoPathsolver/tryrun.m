clear;
global s n m u tpm m0 n1 m1 pm0 pm1 MU vg0 vd0 delta0 v0 alpha0 

%Inputs: s, n, m, u, tpm

tsv0 = 1.0e-2;
%Switching value
ter0 = 1.0e-6;
%Termination tolerance on the value of t
MPmni0 = 2.0e4;
%Maximal number of iterations of the main program
mi0 = 1.0e3;
%Maximal Number of Iterations for corrector step
e0 = 1.0e-1;
de0 = 2;
%Accuracy for starting point of corrector step
hstep0 = 0.1;
%Initial step length
pre0 = 1.0e-15;
pre1 = 1.0e-5;
pre2 = 1.0e-20;
%Control Parameter

delta0 = 0.95;
%Discount factor




s = 1;
n = 2;

m = zeros(s, n);
m(1, 1) = 2;
m(1, 2) = 2;

for q = 1:s
    mp0 = prod(m(q, :));
    pm0(q) = prod(m(q, :));
end

u = zeros(n, mp0, s);
u(:, 1:pm0(1), 1) = [1 0 0 -1; -1 0 0 1];
pv2;

tpm = zeros(n, mp0, s, s);
for q = 1:s
    s1 = m(q, 1);
   
    
    for h = 1:s
        k0 = 0;
        for i1 = 1:s1
           
                for i =1:n
                        tpm(i, k0+1, q, h) = pv0(i1, q, h);
                        k0 = k0 + 1;
                
                end
        end
    end
    
   
end








tic;

n1 = s*n;
m0 = ones(1, s)*m*ones(n, 1);
m1 = m*ones(n, 1);
pm0 = zeros(s, 1);
for q = 1:s
    pm0(q) = prod(m(q, :));
end
pm1 = zeros(s, 1);
for h = 1:s
    pm1(h) = pm0(h)/min(m(h,:));
end

alpha0 = pre0*rand(m0, 1);
%Perturbation vector;

qd0 = m0 + n1 + 1;
%Total Number of Variables

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


%(x,\lambda,mu)

% q0 = 0;
% for q = 1:s
%     k0 = 0;
%     for i = 1:n
%         h0 = q0 + k0;
%         wa = rand(m(q, i), 1);
%         y(h0+1:h0+m(q, i)) = wa/sum(wa);
%         k0 = k0 + m(q, i);
%     end
%     q0 = q0 + m1(q);
% end

%x = pathmcp([y(1:m0)';ones(m0+n1,1)],[zeros(2*m0,1);-inf*ones(n1,1)],inf*ones(2*m0+n1,1),'testexm')
x = pathmcp(rand(m0+n1,1),[zeros(m0,1);-inf*ones(n1,1)],inf*ones(m0+n1,1),'tryexm')
f = tryexm(x,0)

resid = f'*x
time = toc;
fprintf('Computational Time                   = %15.2f\n\n', time);
if abs(resid) > 1.0e-6
  error ('bad solution x from pathmcp run 2');
end

homo;
% x = pathmcp(rand(3,1),zeros(3,1),inf*ones(3,1),'lcp3',zeros(2,3),...
%                    [0; -1],[-1; 1])
% f = lcp3(x,0)
% resid = f'*x
% if abs(resid) > 5e-6
%   error ('bad solution x from pathmcp run 3');
% end

display 'runtests completed OK'

