A = readmatrix("hw0.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
B = readmatrix("B.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
B3 = readmatrix("B3.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
B_subdiag = readmatrix("B_sub.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
C = readmatrix("C.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
C3 = readmatrix("C3.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
C_subdiag = readmatrix("C_sub.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');


B = reshape(B(~isnan(B)), 6, 6);
B3 = B3(~isnan(B3));
B_subdiag = B_subdiag(~isnan(B_subdiag));

C = reshape(C(~isnan(C)), 1200, 1200);
C3 = C3(~isnan(C3));
C_subdiag = C_subdiag(~isnan(C_subdiag));

figure(1);
plot(diag(B), 'X');
title("Diagonal Values of B");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,7]);

figure(2);
plot(-sqrt(diag(B)), 'X');
title("Negative Square Root Diagonal Values of B");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,7]);

figure(3);
plot(B3, 'x');

figure(4);
plot(B_subdiag, 'x');

figure(5);
plot(diag(C), 'X');
title("Diagonal Values of C");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,1201]);

figure(6);
plot(-sqrt(diag(C)), 'X');
title("Negative Square Root Diagonal Values of C");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,1201]);

figure(7);
plot(C3, 'x');

figure(8);
plot(C_subdiag, 'x');