%ASLTP and Square Smoothing

global s n m  m0 m1 vg0 vd0  up0 pp0 v0  pre1 pre2 pre0 de0 e0 hstep0 MPmni0 tsv0 ter0 mi0 beta0

%Inputs: s, n, m, u, tpm

% tsv0 = 1.0e-2;
% %Swirching value
% ter0 = 1.0e-6;
% %Stop criterion on the value of t
 beta0 = 1.0e-12;
% %Smoothing parameter in nu(t)
% MPmni0 = 1.0e4;
% %Maximal number of iterations of the main program
% mi0 = 1.0e3;
% %Maximal Number of Iterations
% e0 = 1.0e-1;
% de0 = 2;
% %Accuracy of starting point in corrector step
% hstep0 = 0.1;
% %Initial step length
% pre0 = 1.0e-15;
% pre1 = 1.0e-5;
% pre2 = 1.0e-20;
% %Control Parameters
% 
% delta0 = 0.95;
% %Discount factor
% 
% %Initialization
% 
% n1 = s*n;
% m0 = ones(1, s)*m*ones(n, 1);
% m1 = m*ones(n, 1);
% pm0 = zeros(s, 1);
% for q = 1:s
%     pm0(q) = prod(m(q, :));
% end
% pm1 = zeros(s, 1);
% for h = 1:s
%     pm1(h) = pm0(h)/min(m(h,:));
% end
% 
% alpha0 = pre0*rand(m0, 1);
% %Perturbation vector;
% 
% qd0 = m0 + n1 + 1;
% %Total Number of Variables
% 
% MU = zeros(m0, n1);
% q0 = 0;
% q1 = 0;
% for q = 1:s
%     k0 = 0;
%     for i = 1:n
%         h0 = q0 + k0;
%         MU(h0+1:h0+m(q, i), q1+i) = 1;
%         k0 = k0 + m(q, i);
%     end
%     q0 = q0 + m1(q);
%     q1 = q1 + n;
% end
% 
% v0 = zeros(m0, 1);
% %Starting point
% q0 = 0;
% for q = 1:s
%     k0 = 0;
%     for i = 1:n
%         h0 = q0 + k0;
%         wb = rand(m(q, i), 1);
%         v0(h0+1:h0+m(q, i)) = wb/sum(wb);
%         k0 = k0 + m(q, i);
%     end
%     q0 = q0 + m1(q);
% end

bv0 = zeros(m0, 1);
%Belief
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        wa = rand(m(q, i), 1);
        bv0(h0+1:h0+m(q, i)) = wa/sum(wa);
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end
up0 = ysguf1(bv0);
pp0 = ysgprob1(bv0);

x = zeros(qd0, 1);
q0 = 0;
q1 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        x(h0+1:h0+m(q, i)) = sqrt(v0(h0+1:h0+m(q, i)));
        x(m0+q1+i) = 0;
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
    q1 = q1 + n;
end
x(qd0) = 2;

dh0 = dltpsgjm1(x);
A1 = dh0*dh0';
R1 = chol(A1);
r = -1;
while r < pre0
    y = rand(qd0, 1);
    w0 = chy(R1, dh0, y);
    vg0 = y - w0;
    r = norm(vg0);
end
vg0 = vg0/norm(vg0);
if vg0(qd0) > 0
    vg0 = -vg0;
end
c0 = sign(det([dh0; vg0']));

time1 = tic;
time0 = cputime;
%Starting time

%Main program starts here
tol0 = 2;
ni = 0;
ni0 = 0;
fcy0 = 1;
while tol0 >= ter0
    tf0 = pre1*tol0;
    %Tolerence of function values
    ops3 = optimset('TolFun', tf0, 'MaxIter', mi0, 'Algorithm', 'trust-region-reflective', 'Display', 'off');
    cer0 = e0*10^(-de0*(1 - tol0/2) - 1/2);
    if (1.005 > tol0) && (tol0 > 0.995)
        time2 = tic;
    end
    if tol0 >= tsv0
        h0 = hstep0;
        %Step length of the predictor step
        y0 = x + h0*vg0;
        rmin = y0(qd0);
        while rmin < pre0
            h0 = 0.618*h0;
            y0 = x + h0*vg0;
            rmin = y0(qd0);
            if h0 < pre2
                rmin = 1;
                fcy0 = -1;
            end
        end
        f0 = dltpsghf1(y0);
        r0 = norm(f0);
        while r0 >= cer0
            h0 = 0.618*h0;
            y0 = x + h0*vg0;
            f0 = dltpsghf1(y0);
            r0 = norm(f0);
            if h0 < pre2
                r0 = -1;
                fcy0 = -1;
            end
        end
        if fcy0 > 0
            vd0 = vg0'*y0;
            z0 = fsolve(@(x) dltpsghf0(x), y0, ops3);
            tol0 = z0(qd0);
            xs = x;
            x = z0;
            dh0 = dltpsgjm1(x);
            A1 = dh0*dh0';
            R1 = chol(A1);
            r = -1;
            while r < pre0
                y = rand(qd0, 1);
                w0 = chy(R1, dh0, y);
                vg0 = y - w0;
                r = norm(vg0);
            end
            vg0 = vg0/norm(vg0);
            c = sign(det([dh0; vg0']));
            if c ~= c0
                vg0 = -vg0;
            end
        else
            tol0 = -1;
        end
    else
        vg0 = (x - xs)/norm(x - xs);
        h0 = hstep0;
        %Step length of the predictor step
        y0 = x + h0*vg0;
        rmin = y0(qd0);
        while rmin < pre0
            h0 = 0.618*h0;
            y0 = x + h0*vg0;
            rmin = y0(qd0);
            if h0 < pre2
                rmin = 1;
                fcy0 = -1;
            end
        end
        f0 = dltpsghf1(y0);
        r0 = norm(f0);
        while r0 >= cer0
            h0 = 0.618*h0;
            y0 = x + h0*vg0;
            f0 = dltpsghf1(y0);
            r0 = norm(f0);
            if h0 < pre2
                r0 = -1;
                fcy0 = -1;
            end
        end
        if fcy0 > 0
            vd0 = vg0'*y0;
            z0 = fsolve(@(x) dltpsghf0(x), y0, ops3);
            tol0 = z0(qd0);
            xs = x;
            x = z0;
        else
            tol0 = -1;
        end
    end
    if tol0 < 1
        ni0 = ni0 + 1;
    end
    ni = ni + 1;
    if ni >= MPmni0
        tol0 = -1;
    end
end

time0 = cputime - time0;
time1 = toc(time1);
time2 = toc(time2);

t0 = x(qd0);
if t0 < 0
    x = xs;
    t0 = x(qd0);
end
y0 = max([x(1:m0)'; zeros(1, m0)])'.^2;

fprintf('___________________________________________________________________________\n');

fprintf('Differentiable ASLTP\n\n');

fprintf('Number of Players                    = %3d\n', n);
fprintf('Number of States                     = %3d\n', s);
fprintf('Number of Pure Strateies of\n');
fprintf('of Player i in State h               =\n');
for q = 1:s
    fprintf('%3d', m(q, :));
    fprintf('\n');
end
fprintf('\n');

fprintf('Number of Iterations                 = %10d\n', ni);
fprintf('Computational Time                   = %15.2f\n', time1);
fprintf('Number of Iterations (t < 1)         = %10d\n', ni0);
fprintf('Computational Time (t < 1)           = %15.2f\n\n', time2);

fprintf('Final value of t                     = %+0.8e\n', t0);
% fprintf('Tolerance for value of t             = %+0.8e\n', ter0);
% fprintf('Switching value                      = %+0.8e\n', tsv0);
% fprintf('Tolerance corrector starting point   = %+0.8e\n', e0);
% fprintf('de0                                  = %3d\n', de0);
% fprintf('Initial Step Length                  = %+0.8e\n', hstep0);
% fprintf('Smoothing Parameter in nu(t)         = %+0.8e\n', beta0);
% fprintf('Computational Time (cpu)             = %15.2f\n', time0);
fprintf('Feasibility                          = %3d\n\n', fcy0);

% fprintf('Stationary Equilibrium');
% q0 = 0;
% for q = 1:s
%     k0 = 0;
%     for i = 1:n
%         h0 = q0 + k0;
%         fprintf('\n');
%         fprintf('x =\n');
%         fprintf('%15.12f %15.12f %15.12f %15.12f %15.12f', y0(h0+1:h0+m(q, i)));
%         k0 = k0 + m(q, i);
%     end
%     q0 = q0 + m1(q);
% end
% 
% fprintf('\n\n');
% 
% fprintf('Payoffs');
% fu = ysguf1(y0);
% fp = ysgprob1(y0);
% q0 = 0;
% for q = 1:s
%     k0 = 0;
%     for i = 1:n
%         h0 = q0 + k0;
%         wp0 = zeros(m(q, i), 1);
%         q1 = 0;
%         for h = 1:s
%             wp0 = wp0 + x(m0+q1+i)*fp(h0+1:h0+m(q, i), h);
%             q1 = q1 + n;
%         end
%         fprintf('\n');
%         fprintf('u = \n');
%         fprintf('%+15.10f %+15.10f %+15.10f %+15.10f %+15.10f', fu(h0+1:h0+m(q, i)) + delta0*wp0);
%         k0 = k0 +  m(q, i);
%     end
%     q0 = q0 + m1(q);
% end
% fprintf('\n');

fprintf('___________________________________________________________________________\n');
