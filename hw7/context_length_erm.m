function correlation = context_length_erm(i, j, l, sigma, nodes, bcs)

if nargin == 5
	delta = (nodes(j, :) - nodes(i, :)).^2;
else
	delta = (nodes(bcs(j), :) - nodes(bcs(i), :)).^2;
end

dist = sqrt(sum(delta));
correlation = sigma^2;

if (l ~= 0)
	correlation = sigma^2 * (1 + dist / l) * exp(-dist / l);
end

end
