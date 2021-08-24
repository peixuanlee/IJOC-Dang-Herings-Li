global s n m  m0  m1  vg0 vd0  v0 alpha0 qd0 n1 timeh

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
alpha0 = pre0*rand(m0, 1);
%Perturbation vector;

qd0 = m0 + n1 + 1;
%Total Number of Variables


% v0 = zeros(m0, 1);
% %Starting point
% q0 = 0;
% for q = 1:s
%     k0 = 0;
%     for i = 1:n
%         h0 = q0 + k0;
%         wa = rand(m(q, i), 1);
%         v0(h0+1:h0+m(q, i)) = wa/sum(wa);
%         k0 = k0 + m(q, i);
%     end
%     q0 = q0 + m1(q);
% end

x = zeros(qd0, 1);
tx0 = sqrt(v0);
x(1:m0) = tx0 - 1;
x(qd0) = 1;
xs = zeros(1,qd0);
for i = 1 : qd0
xs(i) = x(i);
end


dh0 = ycsgjm1(x);
A1 = dh0*dh0';
R1 = chol(A1);
r = -1;
while r <= pre0
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

tic;
%Time counting
%time0 = cputime;

%Main program starts here
tol0 = 1;
ni = 0;
fcy0 = 1;
while tol0 >= ter0
    tf0 = pre1*tol0;
    ops3 = optimset('TolFun', tf0, 'MaxIter', mi0, 'Algorithm', 'trust-region-reflective', 'Display', 'off');
    cer0 = e0*10^(-de0*(1 - tol0));
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
        f0 = ycsghf1(y0);
        r0 = norm(f0);
        while r0 >= cer0
            h0 = 0.618*h0;
            y0 = x + h0*vg0;
            f0 = ycsghf1(y0);
            r0 = norm(f0);
            if h0 < pre2
                r0 = -1;
                fcy0 = -1;
            end
        end
        if fcy0 > 0
            vd0 = vg0'*y0;
            z0 = fsolve(@(x) ycsghf0(x), y0, ops3);
            tol0 = z0(qd0);
            xs = x;
            x = z0;
            dh0 = ycsgjm1(x);
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
        f0 = ycsghf1(y0);
        r0 = norm(f0);
        while r0 >= cer0
            h0 = 0.618*h0;
            y0 = x + h0*vg0;
            f0 = ycsghf1(y0);
            r0 = norm(f0);
            if h0 < pre2
                r0 = -1;
                fcy0 = -1;
            end
        end
        if fcy0 > 0
            vd0 = vg0'*y0;
            z0 = fsolve(@(x) ycsghf0(x), y0, ops3);
            tol0 = z0(qd0);
            xs = x;
            x = z0;
        else
            tol0 = -1;
        end
    end
    ni = ni + 1;
    if ni >= MPmni0
        tol0 = -1;
    end
    
    



end



%time0 = cputime - time0;
timeh = toc;

t0 = x(qd0);
if t0 < 0
    x = xs;
    t0 = x(qd0);
end
tx0 = sqrt(v0)*t0;
w0 = x(1:m0);
r0 = sqrt(w0.^2 + 4*tx0);
y0 = ((w0 + r0)/2).^2;


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
fprintf('Computational Time  of IPM                 = %15.2f\n\n', timeh);


fprintf('Feasibility                          = %3d\n\n', fcy0);


fprintf('Stationary Equilibrium');
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        fprintf('\n');
        fprintf('x =\n');
        fprintf('%15.12f %15.12f %15.12f %15.12f %15.12f', y0(h0+1:h0+m(q, i)));
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end

