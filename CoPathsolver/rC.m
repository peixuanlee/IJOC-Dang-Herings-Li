%More complicated stochastic cases

clear;

nt0 = 20;

global time nn

diary rC.txt;
fprintf('n = %6.2f\n\n', 4);fprintf('s = %6.2f\n\n', 3);fprintf('m = %6.2f\n\n', 5);




time0 = 0; nc = 0;

pd0 = 0.00;
fprintf('Probability density of zeros in payoff matrix = %6.2f\n\n', pd0);

for k = 1:nt0
    fprintf('k = %4d\n', k);
    se425;
    nc = nn + nc;
    time = time + time0;
    time0 = time;
end

fprintf('Successful Times of PathSolver                 = %15.2f\n\n', nc);
fprintf('Average Computational Time of PathSolver                 = %15.2f\n\n', time/nc);

diary off;