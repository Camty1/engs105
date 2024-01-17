function unpacked = packed2unpacked_coeff(packed)
bandwidth = (size(packed,2) - 1) / 2;
num_elems = size(packed,1);
unpacked = zeros(num_elems);
for i=1:2*bandwidth+1
    diag_start = max(1, bandwidth+2-i);
    diag_end = min(num_elems, num_elems + bandwidth+1-i);
    unpacked = unpacked + diag(packed(diag_start:diag_end, i), i-(bandwidth+1));
end
