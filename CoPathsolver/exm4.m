%Example from Sorin (1986)


global s n m u tpm pd0
pd0 = 0.25;

s = 3;
n = 2;

m = zeros(s, n);
m(1, 1) = 2;
m(1, 2) = 2;
m(2, 1) = 1;
m(2, 2) = 1;
m(3, 1) = 1;
m(3, 2) = 1;

u = zeros(n, 4, s);
u(:, :, 1) = [1 0 0 1; 0 2 1 0];
u(:, 1, 2) = [0; 2];
u(:, 1, 3) = [1; 0];

tpm = zeros(2, 4, 3, 3);
tpm(:, :, 1, 1) = [1 1 0 0; 1 0 1 0];
tpm(:, :, 1, 2) = [0 0 1 0; 0 1 0 0];
tpm(:, :, 1, 3) = [0 0 0 1; 0 0 0 1];
tpm(:, 1, 2, 1) = [0; 0];
tpm(:, 1, 2, 2) = [1; 1];
tpm(:, 1, 2, 3) = [0; 0];
tpm(:, 1, 3, 1) = [0; 0];
tpm(:, 1, 3, 2) = [0; 0];
tpm(:, 1, 3, 3) = [1; 1];

trysg;
