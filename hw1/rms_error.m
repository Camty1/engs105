function rms = rms_error(error)
    rms = sqrt(sum(error.^2, "all")/numel(error));
end