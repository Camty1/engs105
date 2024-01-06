%% ENGS 105 HW 0
% Cameron Wolfe 1/9/24

% Reading files
A = readmatrix("hw0.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
B = readmatrix("B.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
B3 = readmatrix("B3.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
B_subdiag = readmatrix("B_sub.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
C = readmatrix("C.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
C3 = readmatrix("C3.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
C_subdiag = readmatrix("C_sub.dat", 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');

% Remove NANs
B = reshape(B(~isnan(B)), 6, 6);
B3 = B3(~isnan(B3));
B_subdiag = B_subdiag(~isnan(B_subdiag));

C = reshape(C(~isnan(C)), 1200, 1200);
C3 = C3(~isnan(C3));
C_subdiag = C_subdiag(~isnan(C_subdiag));

% Plotting settings
set(0,'DefaultFigureWindowStyle','docked');
blue = colororder;
blue = blue(1,:);

% Plotting
figure(1);
plot(diag(B), 'o-', MarkerFaceColor=blue);
title("Diagonal Values of B");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,7]);

figure(2);
plot(-sqrt(diag(B)), 'o-', MarkerFaceColor=blue);
title("Negative Square Root Diagonal Values of B");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,7]);

figure(3);
plot(B3, 'o-', MarkerFaceColor=blue);
title("Third Row of B");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,7]);

figure(4);
plot(B_subdiag, 'o-', MarkerFaceColor=blue);
title("Second Subdiagonal of B");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,5]);

figure(5);
plot(diag(C), 'o-', MarkerFaceColor=blue);
title("Diagonal Values of C");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,1201]);

figure(6);
plot(-sqrt(diag(C)), 'o-', MarkerFaceColor=blue);
title("Negative Square Root Diagonal Values of C");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0,1201]);

figure(7);
plot(C3, 'o-', MarkerFaceColor=blue);
title("Third Row of C");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0, 1201]);

figure(8);
plot(C_subdiag, 'o-', MarkerFaceColor=blue);
title("Second Subdiagonal of C");
xlabel("Index (i)");
ylabel("Magnitude");
xlim([0, 1199]);