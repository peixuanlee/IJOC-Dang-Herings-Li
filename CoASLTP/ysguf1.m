%Computing all u(q, s^i_j, x^{-i}), q = 1,2,...,s, i=1,2,...,n, j = 1,2,...,m(q,i)

function f = ysguf1(x)

global s n m m0 m1 pm0 pm1 u

f = zeros(m0, 1);
q0 = 0;
for q = 1:s
    W0 = zeros(pm1(q), n);
    for i = 1:n
        r0 = 1;
        k0 = 0;
        for j = 1:i-1
            C0 = r0*x(q0+k0+1:q0+k0+m(q, j))';
            r0 = reshape(C0', [], 1);
            k0 = k0 + m(q, j);
        end
        k0 = k0 + m(q, i);
        for j = i+1:n
            C0 = r0*x(q0+k0+1:q0+k0+m(q, j))';
            r0 = reshape(C0', [], 1);
            k0 = k0 + m(q, j);
        end
        W0(1:pm0(q)/m(q, i), i) = r0;
    end
    k0 = 0;
    for i = 1:n
        h0 = pm0(q)/m(q, i);
        f(q0+k0+1:q0+k0+m(q, i)) = reshape(u(i, 1:pm0(q), q), h0, [])'*W0(1:h0, i);
        k0 = k0 + m(q, i);
    end
    q0 = q0 + m1(q);
end
