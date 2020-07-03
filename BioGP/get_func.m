function tree = get_func( Setslog, set_no, index )
%GET_FUNC Summary of this function goes here
%   Detailed explanation goes here
Setslog = importdata(Setslog);
T_set = Setslog.T.set;
digits(3)
path([pwd '\adf'], path);
for i = 1:length(T_set)
    eval(['syms ' T_set{i}]);
end
p = Setslog.dataset(set_no).pareto.P(index);
tree = sym(p.bias);
for i = 1:p.no_roots
    tree = tree + p.root(i).w*eval(cat(2,p.root(i).tree{:}));
end
tree = vpa(tree); tree = expand(tree);
tree = char(tree);
rmpath([pwd '\adf']);
end

