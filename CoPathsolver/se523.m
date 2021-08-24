%clear;
global s n m u tpm pd0 pm0

s = 2;
n = 5;
% pd0 = 0.25;  %desity of utility
% fprintf('Probability density of zeros in payoff matrix = %6.2f\n\n', pd0);


m = zeros(s, n);
for q = 1:s
    for i = 1:n
        m(q, i) = 3;
    end
end

pm0 = zeros(s,1);
for q = 1:s
    mp0 = prod(m(q, :));
    pm0(q) = prod(m(q, :));
end




u = zeros(n, mp0, s);
for q = 1:s
    npf0 = prod(m(q, :));
    u(:, 1:npf0, q)  =  -10 + round(20*rand(n, npf0));
    for i = 1:n
        for j = 1:npf0
            r0 = rand;
            if r0 <= pd0
                u(i, j, q) = 0;
            end
        end
    end
end




pv5;

tpm = zeros(n, mp0, s, s);
for q = 1:s
    s1 = m(q, 1);
    s2 = m(q, 2); 
    s3 = m(q, 3);
    
    for h = 1:s
        k0 = 0;
        for i1 = 1:s1
            for i2 = 1:s2
                 for i3 = 1:s3
                for i =1:n
                        tpm(i, k0+1, q, h) = pv0(i2, i1, q, h);
                        k0 = k0 + 1;
                end
            end
             end
        end
    end
    
   
end
 
trysg;