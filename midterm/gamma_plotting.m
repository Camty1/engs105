N = 101;

gamma_A = @(R) [1 - 3 * R + sqrt(9 * R.^2 - 2 * R + 1);
                1 - 3 * R - sqrt(9 * R.^2 - 2 * R + 1)] / 2;
gamma_B = @(R) [(6 - 4 * R + sqrt(21 * R.^2 - 36 * R - 36)) ./ (5 * R + 12);
                (6 - 4 * R - sqrt(21 * R.^2 - 36 * R - 36)) ./ (5 * R + 12);];

r_values_A = [0.125, 0.25, 0.5, 1];
r_values_B = [0.5, 1, 1.5, 2];

R_values_A = zeros(length(r_values_A), N);
gamma_values_A = zeros(2*length(r_values_A), N);

R_values_B = zeros(length(r_values_B), N);
gamma_values_B = zeros(2*length(r_values_B), N);

figure(1);
tcl = tiledlayout(2,2);
for i=1:length(r_values_A)
    nexttile;
    R_values_A(i,:) = linspace(0, 2*r_values_A(i), N);
    gamma_values_A(2*i-1:2*i, :) = gamma_A(R_values_A(i,:));
    plot(R_values_A(i,:), gamma_values_A(2*i-1:2*i, :));
    hold on;
    plot([0 3*r_values_A(i) nan, 0 3*r_values_A(i)], [-1 -1 nan 1 1], 'k:');
    hold off;
    xlim([0, 2*r_values_A(i)]);
    xlabel("R");
    ylabel("\gamma_A");
    title("Values of \gamma_A for R \in [0, 2r]" + sprintf(", r = %5.3f", r_values_A(i)));
    legend("\gamma^+_A", "\gamma_A^-", "|\gamma_A|=1", "Location", "southwest");
end

figure(2);
tcl = tiledlayout(2,2);
for i=1:length(r_values_B)
    nexttile;
    R_values_B(i,:) = linspace(0, 4*r_values_B(i), N);
    gamma_values_B(2*i-1:2*i, :) = gamma_B(R_values_B(i,:));
    plot(R_values_B(i,:), real(gamma_values_B(2*i-1:2*i, :)));
    hold on;
    plot(R_values_B(i,:), imag(gamma_values_B(2*i-1:2*i, :)));
    plot(R_values_B(i,:), abs(gamma_values_B(2*i-1:2*i, :)));
    plot([0 3*r_values_B(i) nan, 0 3*r_values_B(i)], [-1 -1 nan 1 1], 'k:');
    hold off;
    xlim([0, 2*r_values_B(i)]);
    xlabel("R");
    ylabel("\gamma_B");
    title("Values of \gamma_B for R \in [0, 4r]" + sprintf(", r = %5.3f", r_values_B(i)));
    legend("Re(\gamma_B^+)", "Re(\gamma_B^-)", "Im(\gamma_B^+)", "Im(\gamma_B^-)", "|\gamma_B^+|", "|\gamma_B^-|", "|\gamma_B|=1", "Location", "southwest");
end

figure(3);
tcl = tiledlayout(2,2);
for i=1:length(r_values_B)
    nexttile;
    R_values_B(i,:) = linspace(0, 4*r_values_B(i), N);
    gamma_values_B(2*i-1:2*i, :) = gamma_B(R_values_B(i,:));
    plot(R_values_B(i,:), abs(gamma_values_B(2*i-1:2*i, :)));
    hold on;
    plot([0 5*r_values_B(i) nan, 0 5*r_values_B(i)], [-1 -1 nan 1 1], 'k:');
    hold off;
    xlim([0, 4*r_values_B(i)]);
    xlabel("R");
    ylabel("\gamma_B");
    title("Values of \gamma_B for R \in [0, 4r]" + sprintf(", r = %5.3f", r_values_B(i)));
    legend("|\gamma_B^+|", "|\gamma_B^-|", "|\gamma_B|=1", "Location", "southwest");
end

figure(4);