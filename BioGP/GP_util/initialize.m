function pop = initialize(pop_size, no_roots, spec_depth)

global max_depth Setslog

err_red_lim = Setslog.err_red_lim;

k = 1;
for i = max_depth-spec_depth+1:max_depth
    for j = 1:ceil(pop_size/spec_depth)
        if k <= pop_size
            pop(k).evaluate = 1; pop(k).no_roots = no_roots; pop(k).bias = 0;
            for l = 1:pop(k).no_roots
                if mod(l,2) == mod(k,2)
                    pop(k).root(l) = grow_tree(i, '1');
                else
                    pop(k).root(l) = full_tree(i, '1');
                end
            end
            pop(k).fval = 0; pop(k).cval = 0;
            k = k+1;
        end
    end
end

for i = 1:length(pop)
    for j = 1:pop(i).no_roots
        pop(i).root(j).evaluate = 1; pop(i).root(j).f = []; pop(i).root(j).w = 0;
        tree = pop(i).root(j).tree_index;
        pop(i).root(j).depth = length(dec2ari(max(tree)));
        tree(pop(i).root(j).T.index > 0) = [];
        pop(i).root(j).nodes = length(unique(tree));
        pop(i).root(j).err_red = err_red_lim;
    end
end

pop = pop';

end