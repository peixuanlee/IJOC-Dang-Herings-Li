%Sationary equilibrium system: Differentiable LTP and Square Smoothing

function f = dltpsghf1(x)

global s n n1 m m0 m1 delta0 up0 pp0 beta0 v0 alpha0 

t0 = x(m0+n1+1);

if t0 <= 1
    zeta0 = 0;
else
    zeta0 = 0.25*(2*t0 - 3)^2 + 0.5*(2*t0 - 3) + 0.25;
end

if t0 <= 1 + beta0
    nu0 = t0;
else
    if t0 < 1 + 3*beta0
        nu0 = t0 - (0.25*(t0 - 1 - 2*beta0)^2/beta0 + 0.5*(t0 - 1 - 2*beta0) + 0.25*beta0);
    else
        nu0 = 1 + 2*beta0;
    end
end

y0 = max([x(1:m0)'; zeros(1, m0)])'.^2;
z0 = max([-x(1:m0)'; zeros(1, m0)])'.^2;
u0 = zeros(m0+n1, 1);
u0(1:m0) = y0;
u0(m0+1:m0+n1) = x(m0+1:m0+n1);
fu = ysguf1(u0);
fp = ysgprob1(u0);

f = zeros(m0+n1, 1);
q0 = 0;
q2 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        wp0 = zeros(m(q, i), 1);
        wp1 = zeros(m(q, i), 1);
        q1 = 0;
        for h = 1:s
            wp0 = wp0 + x(m0+q1+i)*fp(h0+1:h0+m(q, i), h);
            wp1 = wp1 + x(m0+q1+i)*pp0(h0+1:h0+m(q, i), h);
            q1 = q1 + n;
        end
        w0 = (1 - nu0)*(fu(h0+1:h0+m(q, i)) + delta0*wp0);
        w1 = nu0*(up0(h0+1:h0+m(q, i)) + delta0*wp1);
        f(h0+1:h0+m(q, i)) = (1 - zeta0)*(w0 + w1) - zeta0*(y0(h0+1:h0+m(q, i)) - v0(h0+1:h0+m(q, i))) + z0(h0+1:h0+m(q, i)) - x(m0+q2+i) - t0*(2 - t0)*alpha0(h0+1:h0+m(q, i));
        f(m0+q2+i) = sum(y0(h0+1:h0+m(q, i))) - 1;
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
    q2 = q2 + n;
end
