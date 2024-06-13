phi = cell(4,1);
for i = 1:4
    phi{i} = @(x) x+i;
end