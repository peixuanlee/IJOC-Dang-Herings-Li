%Example in Herings's paper

global s n m u tpm pm0

s = 2;
n = 2;
m = zeros(2, 2);

m(1, 1) = 2;
m(1, 2) = 2;
m(2, 1) = 2;
m(2, 2) = 2;

pm0 = zeros(s, 1);
for q = 1:s
    pm0(q) = prod(m(q, :));
end

u = zeros(2, 4, 2);
u(:, 1:pm0(1), 1) = [1 0 0 3; -1 0 0 -3];
u(:, 1:pm0(2), 2) = [1 0 0 3; -1 0 0 -3];

tpm = zeros(2, 4, 2, 2);
tpm(:, 1:pm0(1), 1, 1) = [1 0 0 1; 1 0 0 1];
tpm(:, 1:pm0(2), 2, 1) = [1 0 0 1; 1 0 0 1];
tpm(:, 1:pm0(1), 1, 2) = [0 1 1 0; 0 1 1 0];
tpm(:, 1:pm0(2), 2, 2) = [0 1 1 0; 0 1 1 0];

trysg;
