function ind = recat(ind1, ind2, rmv)

global max_arity

ind = ind1;
m = rmv(1)-1; n = rmv(length(rmv))+1;
l = length(ind1.tree_index);

ind.tree = [ind1.tree(1:m) ind2.tree ind1.tree(n:l)];
ind.tree_index = [ind1.tree_index(1:m) ind2.tree_index ind1.tree_index(n:l)];

for i = 1:max_arity
    ind.F(i).index = [ind1.F(i).index(1:m) ind2.F(i).index ind1.F(i).index(n:l)];
end

ind.T.index = [ind1.T.index(1:m) ind2.T.index ind1.T.index(n:l)];

tree = ind.tree_index;
ind.depth = length(dec2ari(max(tree)));
tree(ind.T.index > 0) = [];
ind.nodes = length(unique(tree));
ind.evaluate = 1;

end