%Jacobian matrix

function F = ycsgjm1(x)

global s n m m0 n1 m1 MU delta0 v0 alpha0

t0 = x(m0+n1+1);

tx0 = sqrt(v0)*t0;
dtx0 = sqrt(v0);
w0 = x(1:m0);
r0 = sqrt(w0.^2 + 4*tx0);
y0 = ((r0 + w0)/2).^2;

dy0 = w0 + r0 - 2*tx0./r0;
dz0 = w0 - r0 + 2*tx0./r0;

dy0t = (1 + w0./r0).*dtx0;
dz0t = (1 - w0./r0).*dtx0;

w1 = [y0' x(m0+1:m0+n1)']';
fu = ysguf1(w1);
fp = ysgprob1(w1);
jmu = ysguf2(w1);
jmp = ysgprob2(w1);

pm = zeros(m0, n1);
q0 = 0;
for q = 1:s
    k0 = 0;
    for i = 1:n
        h0 = q0 + k0;
        q1 = 0;
        for h = 1:s
            pm(h0+1:h0+m(q, i), q1+i) = fp(h0+1:h0+m(q, i), h);
            q1 = q1 + n;
        end
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end

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
        z1(h0+1:h0+m(q, i)) = -(fu(h0+1:h0+m(q, i)) + delta0*wp) - 1 - (1 - 2*t0)*alpha0(h0+1:h0+m(q, i));
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end

F11 = (1 - t0)*(jmu + delta0*jmp)*diag(dy0) + diag(dz0);
F12 = -MU + (1 - t0)*delta0*pm;
F13 = z1 + (1 - t0)*(jmu + delta0*jmp)*dy0t + dz0t;

F21 = MU'*diag(dy0);
F22 = zeros(n1);
F23 = MU'*dy0t;

F = [F11 F12 F13; F21 F22 F23];

[U, S, V] =svd(F, 'econ');
F = U*(S + 1.0e-3*t0*eye(m0 + n1))*V';
