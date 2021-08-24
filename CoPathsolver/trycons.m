%Jacobian matrix

function F = trycons(x)

global s n m m0 n1 m1 MU delta0


fp = ysgprob1(x);
jmu = ysguf2(x);
jmp = ysgprob2(x);

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


F11 = jmu + delta0*jmp;
F12 = -MU + delta0*pm;

%F21 = MU';
%F22 = zeros(n1);


%F = [F11 F12; F21 F22];
F = -[F11 F12];


