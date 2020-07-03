function out=f_BioGP(data,Setslog)


plst = {'\adf'};
for i = 1:length(plst)
    path([pwd plst{i}], path);
end
%Setslog=importdata(Setslog);

f_index = Setslog.out_index;
DataSet = Setslog.DataSet;
omin = min(DataSet(:,f_index)); omax = max(DataSet(:,f_index));
in_index = Setslog.in_index;
T_set = Setslog.T.set;

if Setslog.subsets == 1
    p = Setslog.dataset.pareto.P(1);
else
    s = Setslog.subsets + 1;
    p = Setslog.dataset(s).pareto.P(1);
end

input_col=length(Setslog.in_index);
in=data(:,1:input_col);
noexp = length(in(:,1));

for i = 1:length(in_index)
    eval([T_set{i} ' = in(:,i);']);
end
feval = ones(noexp, p.no_roots+1); w = ones(p.no_roots+1,1); w(1) = p.bias;
for i = 1:p.no_roots
    tree = eval(cat(2,p.root(i).tree{:}));
    feval(:,i+1) = tree.*ones(noexp,1);
    w(i+1) = p.root(i).w;
end
traindata2 = feval*w;

out=traindata2;
