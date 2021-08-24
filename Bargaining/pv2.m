%Generating transition probability for 2 players
%pv0(s^2, s^1, q, h) =  trabsition probability from state q to state h 
%when the pure strategy profile (s^2, s^1) is played at state q

rm = 0;
for q = 1:s
    r0 = max(m(q, :));
    if rm < r0
        rm = r0;
    end
end

pv0 = zeros(rm, rm, s, s);
for q = 1:s
    s1 = m(q, 1);
    s2 = m(q, 2);
    for h = 1:s
        pv0(1:s2, 1:s1, q, h) = rand(s2, s1);
    end
    for i1 = 1:s1
        for i2 = 1:s2
                pv0(i2, i1, q, :) = pv0(i2, i1, q, :)/sum(pv0(i2, i1, q, :));
        end
    end
end
