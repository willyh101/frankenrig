function score = avoidHittingCellsToMuch(matrix)

[C,ia,ic] = unique(matrix);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

score = max(value_counts(:,2));
