%Stationary equilibrium system: Logarithmic Barrier (V1) & CY Smoothing

function f = ycsghf1(x)

global s n n1 m m0 m1 delta0 v0 alpha0

t0 = x(m0+n1+1);

tx0 = sqrt(v0)*t0;
w0 = x(1:m0);
r0 = sqrt(w0.^2 + 4*tx0);
y0 = ((w0 + r0)/2).^2;
z0 = ((w0 - r0)/2).^2;

fu = ysguf1(y0);
fp = ysgprob1(y0);

f = zeros(m0+n1, 1);
q0 = 0;
q2 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        wp0 = zeros(m(q, i), 1);
        q1 = 0;
        for h = 1:s
            wp0 = wp0 + x(m0+q1+i)*fp(h0+1:h0+m(q, i), h);
            q1 = q1 + n;
        end
        f(h0+1:h0+m(q, i)) = (1 - t0)*(fu(h0+1:h0+m(q, i)) + delta0*wp0) + z0(h0+1:h0+m(q, i)) - x(m0+q2+i) - t0 - t0*(1-t0)*alpha0(h0+1:h0+m(q, i));
        f(m0+q2+i) = sum(y0(h0+1:h0+m(q, i))) - 1;
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
    q2 = q2 + n;
end
