clear;

nt0 = 10;
%Number of tests

diary A0.txt;

pd0 = 0.00;
fprintf('Probability density of zeros in payoff matrix = %6.2f\n\n', pd0);

for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs225;
end


for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs253;
end


for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs254;
end


for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs255;
end

for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs333;
end

for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs335;
end

for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs425;
end

for k = 1:nt0
    fprintf('k = %4d\n', k);
    inputs525;
end

diary off;

