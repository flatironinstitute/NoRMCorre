function patches = construct_grid_even(grid_size,overlap,dim,min_diff)

N = length(dim);
start_point = cell(N,1);
end_point = cell(N,1);

for i = 1:N
    start_point{i} = [1:grid_size(i):dim(i)-grid_size(i)-overlap(i)+1,max(dim(i)-grid_size(i)-overlap(i)+1,1)];
    end_point{i} = [grid_size(i)+overlap(i):grid_size(i):dim(i),dim(i)];
    if length(start_point{i}) > 1
        if start_point{i}(end) - start_point{i}(end-1) < min_diff(i)        
            start_point{i}(end-1) = [];
            end_point{i}(end-1) = [];
        end
    end
end

start_grid = cell(1,N);
[start_grid{:}] = ndgrid(start_point{:});

end_grid = cell(1,N);
[end_grid{:}] = ndgrid(end_point{:});
patches = zeros(numel(start_grid{1}),2*N);

for i = 1:N
    patches(:,2*i-1) = start_grid{i}(:);
    patches(:,2*i) = end_grid{i}(:);
end