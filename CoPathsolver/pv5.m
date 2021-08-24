%Generating transition probability for 5 players
%pv0(s^5, s^4, s^3, s^2, s^1, q, h) =  trabsition probability from state q to state h 
%when the pure strategy profile (s^5, s^4, s^3, s^2, s^1) is played at stae q
clear sum;
rm = 0;
for q = 1:s
    r0 = max(m(q, :));
    if rm < r0
        rm = r0;
    end
end
pv0 = zeros(rm, rm, rm, rm, rm, s, s);
for q = 1:s
    s1 = m(q, 1);
    s2 = m(q, 2);
    s3 = m(q, 3);
    s4 = m(q, 4);
    s5 = m(q, 5);
    for h = 1:s
        pv0(1:s5, 1:s4, 1:s3, 1:s2, 1:s1, q, h) = rand(s5, s4, s3, s2, s1);
    end
    for i1 = 1:s1
        for i2 = 1:s2
            for i3 = 1:s3
                for i4 = 1:s4
                    for i5 = 1:s5
                        pv0(i5, i4, i3, i2, i1, q, :) = pv0(i5, i4, i3, i2, i1, q, :)/sum(pv0(i5, i4, i3, i2, i1, q, :));
                    end
                end
            end
        end
    end
end
