function [pop F1 F2] = GPevaltree(pop, check_final)

global Setslog setno
global in_index F_bad

chk_w = 1e-10;              %subtrees with absolute weights lower than this are eliminated

for i = 1:length(in_index)
    eval([Setslog.T.set{i} '= Setslog.dataset(setno).in(:,i);']);
end
output = Setslog.dataset(setno).out;
err_red_lim = Setslog.err_red_lim;
lambda = Setslog.lambda;
max_depth = Setslog.max_depth;
if check_final
    for i = 1:length(pop)
        j = 1;
        while j <= pop(i).no_roots
            if pop(i).root(j).err_red < err_red_lim || abs(pop(i).root(j).w) < chk_w
                pop(i).root(j) = [];
                pop(i).no_roots = pop(i).no_roots-1;
                pop(i).evaluate = 1;
            else
                j = j+1;
            end
        end
    end
    j = 1; len = length(pop);
    while j <= len
        if pop(j).no_roots == 0
            pop(j) = []; len = len-1;
        else
            j = j+1;
        end
    end
end
F1 = ones(size(pop)); F2 = F1;
omax = Setslog.dataset(Setslog.no_run).out;
omin = min(omax); omax = max(omax);
for i = 1:length(pop)
    if pop(i).evaluate
        noexp = (length(output));
        
        warning('off', 'all')
        feval = ones(noexp, 1);
        
        for j = 1:pop(i).no_roots
            if pop(i).root(j).err_red < err_red_lim && check_final == 0
                temp_pop = initialize(1, 1, max_depth);
                pop(i).root(j) = temp_pop.root(1);
            end
            if pop(i).root(j).evaluate
                tree = eval(cat(2,pop(i).root(j).tree{:}));
                if isnan(sum(tree)) || isinf(sum(tree))
                    feval(:,j+1) = zeros(noexp,1);
                else
                    feval = [feval tree.*ones(noexp,1)];
                end
                pop(i).root(j).f = feval(:,j+1);
                pop(i).root(j).evaluate = 0;
            else
                feval = [feval pop(i).root(j).f];
            end
        end
        coeff = feval\output;
        
%         feval1 = feval.*(coeff*ones(1,length(feval(:,1))))';
        feval1 = feval(:,2:length(feval(1,:)));
        [Q,~] = qr(feval1,0); auxp = Q\output;
%         [Q,~] = qr(feval,0); auxp = Q\output;
        warning('on', 'all')
        
        depth = ones(1,pop(i).no_roots); nodes = ones(size(depth));
        pop(i).bias = coeff(1);
        for j = 1:pop(i).no_roots
            qi = Q(:,j); auxpp = auxp(j);
%             qi = Q(:,j+1); auxpp = auxp(j+1);
            pop(i).root(j).w = coeff(j+1);
            pop(i).root(j).err_red = (auxpp^2)*(qi'*qi)/(output'*output);
            depth(j) = pop(i).root(j).depth; nodes(j) = pop(i).root(j).nodes;
        end
        out = feval*coeff;
        
        pop(i).fval = sqrt(sum(((out-output)/(omax-omin)).^2)/noexp);
        pop(i).cval = lambda*max(depth) + (1-lambda)*sum(nodes);
        if isnan(pop(i).fval) || isinf(pop(i).fval)
            pop(i).fval = F_bad+eps; pop(i).cval = F_bad+eps;
        end
        % Put constraints on trained output if necesaary
%         if ~isempty(find(out < 0, 1)) || ~isempty(find(out > 1, 1))
%             pop(i).fval = pop(i).fval + 1e-3*(length(find(out < 0)) + length(find(out > 1)));
%         end
        
        pop(i).evaluate = 0;
    end
    F1(i) = pop(i).fval; F2(i) = pop(i).cval;
end

end