function pop = mutate(pop)

global F T
global max_arity
global Pm_stand Pm_small Pm_mono

for i = 1:length(pop)
    ind = pop(i);
    flag = 1; r = rand;
    while flag
        if r <= Pm_stand
            %standard mutation
            ind.evaluate = 1;
            root_id = ceil(rand*ind.no_roots); ind_new = ind.root(root_id); ind_new.evaluate = 1;
            r = ind_new.tree_index;
            d = ceil(rand*ind_new.depth);
            r = unique(r(r >= max_arity^(d-1) & r <= max_arity^d-1));
            r = r(ceil(rand*length(r)));
            index = dec2ari(r);
            root_new = grow_tree(d, index);
            s = find(ind_new.tree_index == r);
            ind.root(root_id) = recat(ind_new, root_new, s); flag = 0;
        elseif r <= Pm_small + Pm_stand
            %small mutation
            ind.evaluate = 1;
            root_id = ceil(rand*ind.no_roots); ind_new = ind.root(root_id); ind_new.evaluate = 1;
            r = ind_new.tree_index;
            d = ceil(rand*ind_new.depth);
            r = unique(r(r >= max_arity^(d-1) & r <= max_arity^d-1));
            r = r(ceil(rand*length(r)));
            d = find(ind_new.tree_index == r);
            arity = length(d) - 1; d = d(1);
            if arity > 0
                r = find(F.sets == arity);
                if length(r) > 1
                    r(strcmp(F.set(r), ind_new.tree{d}(1:length(ind_new.tree{d})-1))) = [];
                end
                r = r(ceil(rand*length(r)));
                ind_new.tree(d) = {[F.set{r} '(']};
            else
                r = ind_new.tree{d};
                if isnan(str2double(r))
                    r = 1:length(T.set);
                    r(strcmp(T.set, 'ERC')) = []; r(strcmp(T.set, 'erc')) = [];
                    if length(r) > 1
                        r(strcmp(T.set(r), ind_new.tree{d})) = [];
                    end
                    r = r(ceil(rand*length(r)));
                    ind_new.tree(d) = T.set(r);
                else
                    ind_new.tree(d) = {num2str(str2double(r)+randn)};
                end
            end
            ind.root(root_id) = ind_new;  flag = 0;
        elseif r <= Pm_mono + Pm_small + Pm_stand
            %monoparental xover
            ind.evaluate = 1;
            root = struct2cell(ind.root(:)); names = fieldnames(ind.root);
            index = 1:ind.no_roots; depth = cell2mat(root(strcmp(names,'depth'),:)); index(depth == 1) = [];
            if isempty(index)
                r = rand*(Pm_small + Pm_stand); continue
            end
            root_id = index(ceil(rand(1,2)*length(index))); par = ind.root(root_id);
            if root_id(1) ~= root_id(2)
                par = xover(par, 1);
            else
                par = xover(par, 2); root_id = root_id(1);
            end
            ind.root(root_id) = par;
            flag = 0;
        else
            flag = 0;
        end
    end
    pop(i) = ind;
end

end