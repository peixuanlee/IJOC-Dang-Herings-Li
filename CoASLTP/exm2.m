%Example from Thuijsman

global s n m u tpm

s = 4;
n = 2;

m = zeros(s, n);
m(1, 1) = 2;
m(1, 2) = 2;
m(2, 1) = 2;
m(2, 2) = 2;
m(3, 1) = 1;
m(3, 2) = 1;
m(4, 1) = 1;
m(4, 2) = 1;
% 
 u = zeros(n, 4, s);
% u(:, :, 1) = [0 0 0 0; 0 0 0 0];
% u(:, :, 2) = [0 0 0 0; 0 0 0 0];
% u(:, 1, 3) = [1; -1];
% u(:, 1, 4) = [-1; 1];
% 
% tpm = zeros(2, 4, 4, 4);
% tpm(:, :, 1, 1) = [1 0 0 0; 1 0 0 0];
% tpm(:, :, 1, 2) = [0 0 0 1; 0 0 0 1];
% tpm(:, :, 1, 3) = [0 1 1 0; 0 1 1 0];
% tpm(:, :, 1, 4) = [0 0 0 0; 0 0 0 0];
% tpm(:, :, 2, 1) = [1 0 0 0; 1 0 0 0];
% tpm(:, :, 2, 2) = [0 0 0 1; 0 0 0 1];
% tpm(:, :, 2, 3) = [0 0 0 0; 0 0 0 0];
% tpm(:, :, 2, 4) = [0 1 1 0; 0 1 1 0];
% tpm(:, 1, 3, 1) = [0; 0];
% tpm(:, 1, 3, 2) = [0; 0];
% tpm(:, 1, 3, 3) = [1; 1];
% tpm(:, 1, 3, 4) = [0; 0];
% tpm(:, 1, 4, 1) = [0; 0];
% tpm(:, 1, 4, 2) = [0; 0];
% tpm(:, 1, 4, 3) = [0; 0];
% tpm(:, 1, 4, 4) = [1; 1];

u(:, :, 1) = [1 0 -7 1; 1 0 -7 1];
u(:, :, 2) = [1 0 0 3; -1 0 0 -3];
u(:, 1, 3) = [0; 0];
u(:, 1, 4) = [0; 0];

tpm = zeros(2, 4, 4, 4);
tpm(:, :, 1, 1) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, :, 1, 2) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, :, 1, 3) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, :, 1, 4) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, :, 2, 1) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, :, 2, 2) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, :, 2, 3) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, :, 2, 4) = [1 1 1 1; 1 1 1 1]/4;
tpm(:, 1, 3, 1) = [1; 1]/4;
tpm(:, 1, 3, 2) = [1; 1]/4;
tpm(:, 1, 3, 3) = [1; 1]/4;
tpm(:, 1, 3, 4) = [1; 1]/4;
tpm(:, 1, 4, 1) = [1; 1]/4;
tpm(:, 1, 4, 2) = [1; 1]/4;
tpm(:, 1, 4, 3) = [1; 1]/4;
tpm(:, 1, 4, 4) = [1; 1]/4;
tpm(:, 1, 4, 4) = [1; 1]/4;



ycsgse;
