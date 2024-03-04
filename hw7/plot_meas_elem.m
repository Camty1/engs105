elem = readmatrix("hw44.ele", "FileType", "text");
elem = elem(:, 2:4);
node = readmatrix("hw44.nod", "FileType", "text");
node = node(:, 2:3);
sample = readmatrix("sample_points.dat");
program_elems = readmatrix("output/meas_elems.dat");

positions = zeros(2,4);

co = colororder;

close all;

for i=1:12
	for j=1:3
		positions(:,j) = node(elem(program_elems(i), j), :)';
	end
	positions(:, 4) = positions(:, 1);
	co_index = mod(i-1, length(co)) + 1;
	plot(sample(i, 1), sample(i, 2), '.', "Color", co(co_index, :), "DisplayName", ['Measurement: ', num2str(i)]);
	if i == 1
		hold on;
	end
	plot(positions(1,:)', positions(2,:)', "Color", co(co_index, :), "DisplayName", ['Element: ', num2str(program_elems(i))]);
	
end

title("Measurements and Elements");
xlabel("x");
ylabel("y");
legend("Location", "eastoutside");

hold off;
