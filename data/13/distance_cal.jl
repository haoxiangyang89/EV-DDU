using CSV, DataFrames;
# Calculate the distance between buses and demand sectors
coordinate_list_i = [(-600, 0),
                     (0, 2000),
                     (-500, 2000),
                     (-800, 2000),
                     (500, 2000),
                     (500, 2000),
                     (0, 4000),
                     (-300, -800),
                     (0, 0),
                     (500, 0),
                     (0, -1000),
                     (-300, 0),
                     (0, 0)];
coordinate_list_j = [(-300, 3000),
                     (400, 2500),
                     (-500, 1000),
                     (-200, -800),
                     (400, -300),
                     (0, -1100),
                     (-1000, 200),
                     (800,3800),
                     (700,400),
                     (500,1500)];
distance_ij = zeros(length(coordinate_list_i), length(coordinate_list_j));
for i in eachindex(coordinate_list_i)
    for j in eachindex(coordinate_list_j)
        distance_ij[i, j] = sqrt((coordinate_list_i[i][1] - coordinate_list_j[j][1])^2 +
                                 (coordinate_list_i[i][2] - coordinate_list_j[j][2])^2)
    end
end
df = DataFrame(distance_ij, :auto);
CSV.write("./data/13/case13_rdis_new.csv",df);