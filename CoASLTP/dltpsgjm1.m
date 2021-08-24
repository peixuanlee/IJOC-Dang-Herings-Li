%Jacobian matrix

function F = dltpsgjm1(x)

global s n m m0 n1 m1 MU delta0 up0 pp0 v0 beta0 alpha0

t0 = x(m0+n1+1);

if t0 <= 1
    zeta0 = 0;
    dzeta0 = 0;
else
    zeta0 = 0.25*(2*t0 - 3)^2 + 0.5*(2*t0 - 3) + 0.25;
    dzeta0 = 2*(t0 - 1);
end
if t0 <= 1 + beta0
    nu0 = t0;
    dnu0 = 1;
else
    if t0 < 1 + 3*beta0
        nu0 = t0 - (0.25*(t0 - 1 - 2*beta0)^2/beta0 + 0.5*(t0 - 1 - 2*beta0) + 0.25*beta0);
        dnu0 = 1 - (0.5*(t0 - 1 - 2*beta0)/beta0 + 0.5);
    else
        nu0 = 1 + 2*beta0;
        dnu0 = 0;
    end
end

y0 = max([x(1:m0)'; zeros(1, m0)])';
z0 = min([x(1:m0)'; zeros(1, m0)])';

u0 = zeros(m0+n1, 1);
u0(1:m0) = y0.^2;
u0(m0+1:m0+n1) = x(m0+1:m0+n1);
fu = ysguf1(u0);
jmu = ysguf2(u0);
fp = ysgprob1(u0);
jmp = ysgprob2(u0);

f0 = zeros(m0, 1);
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        wp1 = zeros(m(q, i), 1);
        q1 = 0;
        for h = 1:s
            wp1 = wp1 + x(m0+q1+i)*pp0(h0+1:h0+m(q, i), h);
            q1 = q1 + n;
        end
        f0(h0+1:h0+m(q, i)) = up0(h0+1:h0+m(q, i)) + delta0*wp1;
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end

pm = zeros(m0, n1);
pm0 = zeros(m0, n1);
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        q1 = 0;
        for h = 1:s
            pm(h0+1:h0+m(q, i), q1+i) = fp(h0+1:h0+m(q, i), h);
            pm0(h0+1:h0+m(q,i), q1+i) = pp0(h0+1:h0+m(q, i), h);
            q1 = q1 + n;
        end
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end
pm = delta0*pm;
pm0 = delta0*pm0;

z1 = zeros(m0, 1);
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        wx = zeros(s, 1);
        qx = 0;
        for h = 1:s
            wx(h) = x(m0+qx+i);
            qx = qx + n;
        end
        wp = fp(h0+1:h0+m(q, i), :)*wx;
        wp0 = pp0(h0+1:h0+m(q, i), :)*wx;
        w0 = fu(h0+1:h0+m(q, i)) + delta0*wp;
        w1 = up0(h0+1:h0+m(q, i)) + delta0*wp0;
        w2 = u0(h0+1:h0+m(q, i)) - v0(h0+1:h0+m(q, i));    
        z1(h0+1:h0+m(q, i)) = -dzeta0*((1 - nu0)*w0 + nu0*w1) + (1 - zeta0)*(-dnu0*w0 + dnu0*w1) - dzeta0*w2 - 2*(1 - t0)*alpha0(h0+1:h0+m(q, i));        
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end

F1 = [((1 - zeta0)*((1 - nu0)*(jmu + delta0*jmp)*diag(2*y0)) - zeta0*diag(2*y0) + diag(2*z0)) (-MU + (1 - zeta0)*((1 - nu0)*pm + nu0*pm0)) z1];
F2 = [MU'*diag(2*y0) zeros(n1, n1+1)];
F = [F1; F2];

[U, S, V] =svd(F, 'econ');
F = U*(S + 1.0e-3*t0*eye(m0 + n1))*V';
