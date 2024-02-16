lhs_band = readmatrix("output/banded_coeff.dat");
rhs = readmatrix("output/rhs.dat");

lhs = zeros(size(rhs));
bandwidth = (size(lhs_band, 2) - 1) / 2;

for i=0:bandwidth
    if i == 0
        lhs = lhs + diag(lhs_band(:, bandwidth + 1));
    else
        lhs = lhs + diag(lhs_band(1:end-i, bandwidth + 1 + i), i);
        lhs = lhs + diag(lhs_band(i+1:end, bandwidth + 1 - i), -i);
    end

end

u = lhs \ rhs