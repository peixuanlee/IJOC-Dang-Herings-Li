%Find a predictor direction

function w0 = chy(R1, dh, y)

w1 = dh*y;
y1 = linsolve(R1', w1);
y2 = linsolve(R1, y1);
w0 = dh'*y2;