function ind = full_tree(depth, index)
%full_tree modified
global F T
global max_depth max_arity

d = ari2dec(index);
if depth >= max_depth
    i = ceil(rand*length(T.set));
    if strcmp(T.set{i}, 'ERC') || strcmp(T.set{i}, 'erc')
        ind.tree = {num2str(randn)};
    else
        ind.tree = T.set(i);
    end
    ind.tree_index = d;
    ind.T.index = d;
    f = 0;
else
    i = ceil(rand*length(F.sets));
    f = F.sets(i);
    ind.tree = {[F.set{i} '(']};
    ind.tree_index = d;
    ind.T.index = 0;
    temp = 0;
    for i = 1:f
        child(i) = full_tree(depth+1, [index dec2ari(temp)]);
        ind.tree = [ind.tree child(i).tree {','}];
        ind.tree_index = [ind.tree_index child(i).tree_index d];
        ind.T.index = [ind.T.index child(i).T.index 0];
        temp = temp + 1;
    end
    ind.tree(length(ind.tree)) = {')'};
end

for j = 1:max_arity
    if j == f
        i = d;
    else
        i = 0;
    end
    temp = i;
    for k = 1:f
        temp = [temp child(k).F(j).index i];
    end
    ind.F(j).index = temp;
end

end