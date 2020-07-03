function ind = rename_ind(ind, depth, new_root)

global max_arity

for i = 1:length(ind.tree_index)
    t = dec2ari(ind.tree_index(i));
    t(1:depth) = [];
    ind.tree_index(i) = ari2dec([new_root t]);
end
for i = 1:max_arity
    t = find(ind.F(i).index);
    ind.F(i).index(t) = ind.tree_index(t);
end

t = find(ind.T.index);
ind.T.index(t) = ind.tree_index(t);
ind.evaluate = 1;

end

