%five fundamental examples

clear;

nt0 = 100;

global timep nn timeh pd0

diary com.txt;

time0 = 0; nc = 0; time1 = 0;

pd0 = 0.5;
fprintf('Probability density of zeros in payoff matrix = %6.2f\n\n', pd0);

for k = 1:nt0
    fprintf('k = %4d\n', k);
    se225;
    nc = nn + nc;
    timep = timep + time0;
    time0 = timep;
    timeh = timeh +time1;
    time1 = timeh;
end

fprintf('Successful Times of PathSolver                 = %15.2f\n\n', nc);
fprintf('Average Computational Time of PathSolver                 = %15.2f\n\n', timep/nc);
%fprintf('Average Computational Time of IPM                 = %15.2f\n\n', timeh/nt0);

diary off;