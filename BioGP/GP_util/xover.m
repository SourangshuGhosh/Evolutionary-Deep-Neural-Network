function pop_new = xover(pop, mono)

global max_arity max_depth
global choose_xfunc Px_stand Px_fair

dox = zeros(floor(length(pop(:,1))/2), 2);
pop_new = pop;

if mono == 0
    if length(pop) == 2
        dox = [1 2];
    else
        index = 1:length(pop(:,1)); k = 1;
        while ~isempty(index)
            if length(index) >= 2
                for j = 1:2
                    select_index = ceil(rand*length(index));
                    dox(k,j) = index(select_index);
                    index(select_index) = [];
                end
            else
                pop_new(2*k-1) = pop(index);
                index = [];
            end
            k = k + 1;
        end
    end
    
    for i = 1:length(dox(:,1))
        
        pop_new(2*i-1) = pop(dox(i,1)); pop_new(2*i) = pop(dox(i,2));
        root_id = ceil(rand(1,2).*[pop_new(2*i-1).no_roots pop_new(2*i).no_roots]);
        p = [pop_new(2*i-1).root(root_id(1)); pop_new(2*i).root(root_id(2))];
        
        depth = [p(1).depth p(2).depth]; node = ones(1,2); r = rand;
        if r <= Px_fair
            %Height fair xover
            d = ones(1,2)*ceil(rand*min(depth));
            x_type = 1;
            pop_new(2*i-1).evaluate = 1;
            pop_new(2*i).evaluate = 1;
        elseif r <= Px_fair + Px_stand
            %Koza Standard Xover
            d(1) = ceil(rand*depth(1)); r = 1:depth(2); r(r==d(1)) = [];
            if ~isempty(r)
                d(2) = r(ceil(rand*length(r))); x_type = 2;
                if d(2) < d(1), k = [2 1]; d = d(k); p = p(k); depth = depth(k); end
            else
                d(2) = d(1); x_type = 1;
            end
            pop_new(2*i-1).evaluate = 1;
            pop_new(2*i).evaluate = 1;
        else
            continue
        end
        
        k = 1;
        while k <= 2
            r = rand;
            if r < choose_xfunc && d(k) ~= depth(k)
                temp = struct2cell(p(k).F(:)); names = fieldnames(p(k).F);
                temp = cell2mat(temp(strcmp(names, 'index'),:));
            else
                temp = p(k).T.index;
            end
            r = temp;
            temp = unique(r(r >= max_arity^(d(k)-1) & r <= max_arity^d(k)-1));
            if isempty(temp)
                temp = struct2cell(p(k).F(:)); names = fieldnames(p(k).F);
                temp = cell2mat(temp(strcmp(names, 'index'),:));
                r = temp;
                temp = unique(r(r >= max_arity^(d(k)-1) & r <= max_arity^d(k)-1));
            end
            node(k) = temp(ceil(rand*length(temp)));
            s = find(p(k).tree_index == node(k));
            depn(k) = length(dec2ari(max(p(k).tree_index(s(1):s(length(s))))));
            if k == 2 && x_type == 2
                if (max_depth-d(2)) < (depn(1)-d(1)) || (max_depth-d(1)) < (depn(2)-d(2))
                    d(2) = d(2) - 1; k = 1; if d(1) == d(2), x_type = 1; end
                end
            end
            k = k+1;
        end
        
        offspring = exchange(p, d, node);
        pop_new(2*i-1).root(root_id(1)) = offspring(1);
        pop_new(2*i).root(root_id(2)) = offspring(2);
    end
    
else
    p = pop; depth = [p(1).depth p(2).depth]; node = zeros(1,2);
    if mono == 2
        d = ceil(ceil(rand(1,2).*(depth-2)+1)+eps);
    else
        d = ceil(rand(1,2).*(depth-1)+1);
    end
    if d(2) < d(1), k = [2 1]; d = d(k); p = p(k); depth = depth(k); end 
    
    k = 1;
    while k <= 2
        if d(k) ~= depth(k)
            temp = struct2cell(p(k).F(:)); names = fieldnames(p(k).F);
            temp = cell2mat(temp(strcmp(names, 'index'),:));
        else
            temp = p(k).T.index;
        end
        r = temp; temp = unique(r(r >= max_arity^(d(k)-1) & r <= max_arity^d(k)-1));
        node(k) = temp(ceil(rand*length(temp))); node_ari(k) = {dec2ari(node(k))};
        s = find(p(k).tree_index == node(k));
        depn(k) = length(dec2ari(max(p(k).tree_index(s(1):s(length(s))))));
        if k == 2
            if mono == 2 && length(node_ari{1}) < length(node_ari{2}) && strcmp(node_ari{1}, node_ari{2}(1:length(node_ari{1})))
                temp = p(2).tree_index;
                s = find(temp == node(2)); temp(s(1):s(length(s))) = []; s = find(temp == node(1));
                s = length(dec2ari(max(temp(s(1):s(length(s))))));
                if s + 1 > max_depth, d(2) = d(2) - 1; k = 1; end
            elseif ((max_depth-d(2)) < (depn(1)-d(1)) || (max_depth-d(1)) < (depn(2)-d(2)))
                d(2) = d(2) - 1; k = 1;
            end
        end
        k = k+1;
    end
    
    if mono == 2
        if node(1) < node(2)
            k = [2 1]; d = d(k); node = node(k); p = p(k); node_ari = node_ari(k);
        end
        if node(1) == node(2)
            offspring = p;
        elseif strcmp(node_ari{2}, node_ari{1}(1:length(node_ari{2})))
            offspring = exchange(p, d, node); r = find(offspring(2).tree_index == node(2));
            r = offspring(2).tree_index(r(1)+1:r(length(r))-1); r = r(r>= max_arity^d(2) & r<= max_arity^(d(2)+1)-1);
            r = unique(r); node_new = r(ceil(rand*length(r)));
            offspring = exchange([offspring(2) p(1)], [d(2)+1 d(1)], [node_new node(1)]);
            offspring = exchange(offspring, [d(2)+1 d(2)], [node_new node(2)]);
        else
            offspring = exchange(p, d, node);
            offspring = exchange(offspring, d([2 1]), node([2 1]));
        end
        offspring = offspring(1);
    else
        offspring = exchange(p, d, node);
    end
    pop_new = offspring;
end

end

function off = exchange(par, d, node)
j = 1; off = par;
for k = 1:2
    r = {d(k) dec2ari(node(k+j))};
    s = find(par(k).tree_index == node(k));
    sub(k) = snatch_tree(par(k), s);
    if node(1) ~= node(2)
        sub(k) = rename_ind(sub(k), r{1}, r{2});
    end
    s = find(par(k+j).tree_index == node(k+j));
    off(k+j) = recat(par(k+j), sub(k), s);
    j = -1;
end
end