function[x1 x2 neighbors] = cellHelper(x1 x2 neighbors);
    %Takes in x1 x2 and neighbors in format of allWithPairs function and creates cells that contain all x1 x2 and neighbor info, not just for particles with index higher.
    n = size(neighbors, 1);
    for i = 1:n
        m = size(neighbors{i},2);
        for j = 1:m
            neighbors{neighbors{i}(j)} = [neighbors{neighbors{i}(j)}(1:i - 1) i neighbors{neighbors{i}(j)}(i + 1:end)];

            x1{neighbors{i}(j)} = [x1{neighbors{i}(j)}(1:i - 1)  x1{neighbors{i}(j)}(i + 1:end)];

            x2{neighbors{i}(j)} = [x2{neighbors{i}(j)}(1:i - 1) i x2{neighbors{i}(j)}(i + 1:end)];
            

end