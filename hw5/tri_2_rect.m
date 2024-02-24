element_list = readmatrix("epeltr4.dat");
node_list = readmatrix("npeltr4.dat");
heat_rate_list = readmatrix("bpelt4.dat", "Delimiter", ",");

disp(heat_rate_list(1:110,:));

rect_count = 0;
tri_count = 0;

tri_element_list = [];
tri_heat_rate_list = [];

for i=1:length(element_list)
    if (element_list(i, 4) ~= element_list(i, 5))
        rect_count = rect_count + 1;
        tri_element_list = [tri_element_list; element_list(i, [1 2 3 4 4 6]); element_list(i, [1 2 4 5 5 6])];
        tri_heat_rate_list = [tri_heat_rate_list; heat_rate_list(i, :); heat_rate_list(i, :)];
    else
        tri_count = tri_count + 1;
        tri_element_list = [tri_element_list; element_list(i,:)];
        tri_heat_rate_list = [tri_heat_rate_list; heat_rate_list(i,:)];
    end
end

for i=1:size(tri_element_list, 1)
    tri_element_list(i, 1) = i;
    tri_heat_rate_list(i, 1) = i;
    tri_heat_rate_list(i, 2) = i;
end

disp(rect_count);
disp(tri_count);
disp(2 * rect_count + tri_count);
disp(size(tri_element_list));
disp(tri_heat_rate_list(1:100, :));

writematrix(tri_element_list, "epeltr4_tri.dat");
writematrix(tri_heat_rate_list, "bpelt4_tri.dat");
