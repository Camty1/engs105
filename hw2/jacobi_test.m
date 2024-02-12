A = [3 1; 2 4];
b = [-1; 1];
x_true = A\b;
x = zeros(2,1);

err = max(abs(x-x_true));
A_diag = diag(diag(A));
below = [0 0; 2 0];
above = [0 1; 0 0];
while err > 1e-5
    disp(err);
    x = A_diag\(-(above+below)*x + b);
    err = max(abs(x-x_true));
end