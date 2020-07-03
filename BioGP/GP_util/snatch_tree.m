function ind = snatch_tree(ind, index)

global max_arity

m = index(1); n = index(length(index));
ind.tree = ind.tree(m:n);
ind.tree_index = ind.tree_index(m:n);
for i = 1:max_arity
    ind.F(i).index = ind.F(i).index(m:n);
end

ind.T.index = ind.T.index(m:n);
ind.evaluate = 1;

end