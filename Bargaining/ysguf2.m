%Jacoban matrix of (u(q, s^i_j, x^{-i}): q = 1,2,...,s, i = 1,2,...,n, j = 1,2,...,m(q,i))

function fm = ysguf2(x)

global s n m m0 m1 pm0 u

fm = zeros(m0, m0);

q0 = 0;
for q = 1:s
    h4 = 0;
    for i = 1:n
        h0 = pm0(q)/m(q, i);
        h5 = 0;
        for j = 1:i-1
            k0 = 0;
            r0 = 1;
            for h = 1:j-1
                C0 = r0*x(q0+k0+1:q0+k0+m(q, h))';
                r0 = reshape(C0', [], 1);
                k0 = k0 + m(q, h);
            end
            k0 = k0 + m(q, j);
            for h = j+1:i-1
                C0 = r0*x(q0+k0+1:q0+k0+m(q, h))';
                r0 = reshape(C0', [], 1);
                k0 = k0 + m(q, h);
            end
            k0 = k0 + m(q, i);
            for h = i+1:n
                C0 = r0*x(q0+k0+1:q0+k0+m(q, h))';
                r0 = reshape(C0', [], 1);
                k0 = k0 + m(q, h);
            end
            u1 = reshape(u(i, 1:pm0(q), q), h0, [])';
            h1 = 1;
            if i-1 >= j+1
                h1 = prod(m(q, j+1:i-1));
                if i+1 <= n
                    h1 = h1*prod(m(q, i+1:n));
                end
            else
                if i+1 <= n
                    h1 = prod(m(q, i+1:n));
                end
            end
            for k = 1:m(q, i)
                u2 = reshape(u1(k, :), h1, []);
                h2 = h0/h1;
                h3 = h2/m(q, j);
                u3 = zeros(h1, h3);
                u4 = zeros(h1*h3, m(q, j));
                for ka = 1:m(q, j)
                    for kb = 1:h3
                        u3(:, kb) = u2(:, (kb-1)*m(q, j)+ka);
                    end
                    u4(:, ka) = reshape(u3, [], 1);
                end
                fm(q0+h4+k, q0+h5+1:q0+h5+m(q, j))= r0'*u4;
            end
            h5 = h5 + m(q, j);
        end
        h5 = h5 + m(q, i);
        for j = i+1:n
            k0 = 0;
            r0 = 1;
            for h = 1:i-1
                C0 = r0*x(q0+k0+1:q0+k0+m(q, h))';
                r0 = reshape(C0', [], 1);
                k0 = k0 + m(q, h);
            end
            k0 = k0 + m(q, i);
            for h = i+1:j-1
                C0 = r0*x(q0+k0+1:q0+k0+m(q, h))';
                r0 = reshape(C0', [], 1);
                k0 = k0 + m(q, h);
            end
            k0 = k0 + m(q, j);
            for h = j+1:n
                C0 = r0*x(q0+k0+1:q0+k0+m(q, h))';
                r0 = reshape(C0', [], 1);
                k0 = k0 + m(q, h);
            end
            u1 = reshape(u(i, 1:pm0(q), q), h0, [])';
            h1 = 1;
            if j+1 <= n
                h1 = prod(m(q, j+1:n));
            end
            for k = 1:m(q, i)
                u2 = reshape(u1(k, :), h1, []);
                h2 = h0/h1;
                h3 = h2/m(q, j);
                u3 = zeros(h1, h3);
                u4 = zeros(h1*h3, m(q, j));
                for ka = 1:m(q, j)
                    for kb = 1:h3
                        u3(:, kb) = u2(:, (kb-1)*m(q, j)+ka);
                    end
                    u4(:, ka) = reshape(u3, [], 1);
                end
                fm(q0+h4+k, q0+h5+1:q0+h5+m(q, j))= r0'*u4;
            end
            h5 = h5 + m(q, j);
        end
        h4 = h4 + m(q, i);
    end
    q0 = q0 + m1(q);
end
    