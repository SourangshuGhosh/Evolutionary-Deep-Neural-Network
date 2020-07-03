function mico_table = get_trained_data(Setslog, set_no, yl,savedir)

global ploton saveon fig

plst = {'\adf'};
for i = 1:length(plst)
    path([pwd plst{i}], path);
end

f_index = Setslog.out_index;
DataSet = Setslog.DataSet;
omin = min(DataSet(:,f_index)); omax = max(DataSet(:,f_index));
in_index = Setslog.in_index;
T_set = Setslog.T.set;

p = Setslog.dataset(set_no).pareto.P(1);

mico_table = zeros(1,Setslog.no_run);

for l = 1:Setslog.no_run
    in = Setslog.dataset(l).in;
    out = Setslog.dataset(l).out;
    noexp = length(in(:,1));

    for i = 1:length(in_index)
        eval([T_set{i} ' = in(:,i);']);
    end
    feval = ones(noexp, p.no_roots+1); w = ones(p.no_roots+1,1); w(1) = p.bias;
    for i = 1:p.no_roots
        tree = eval(cat(2,p.root(i).tree{:}));
        %         if isnan(sum(tree)) || isinf(sum(tree))
        %             feval(:,i+1) = zeros(noexp,1);
        %         else
        feval(:,i+1) = tree.*ones(noexp,1);
        w(i+1) = p.root(i).w;
        %     end
    end
    traindata2 = feval*w;
    
    mico_table(1,l) = sqrt(sum(((out-traindata2)/(omax-omin)).^2)/noexp);
    
    traindata1 = Setslog.dataset(l).data_index;
    if l == set_no && l == Setslog.no_run, ploton = 1; end
    if ploton
        set(0,'CurrentFigure',fig(1)); clf
        plot(traindata1, DataSet(traindata1,f_index), '-rd'); hold on
        if l == set_no
            plot(traindata1, traindata2, '--ks')
            legend('Experimental Data', ['Trained Data by GPtree ' num2str(l)])
        else
            plot(traindata1, traindata2, '--cs')
            legend('Experimental Data', ['Predicted Data by GPtree ' num2str(set_no) ' on dataset ' num2str(l)])
        end
        xlabel('data index')
        ylabel(yl)
        pause(1)
        set(0,'CurrentFigure',fig(2)); clf
        plot(DataSet(traindata1,f_index), DataSet(traindata1,f_index), '-kd', 'LineWidth', 2); hold on
        pol = polyfit(DataSet(traindata1,f_index), traindata2, 1);
        text(1.5*0.5*(omin+omax), 0.5*0.5*(omin+omax), ['slope of fitted line: ' num2str(pol(1))], 'Fontsize', 12)
        pol = pol(1)*DataSet(traindata1,f_index) + pol(2);
        plot(DataSet(traindata1,f_index), pol, '.-b', 'LineWidth', 2);
        xlabel(['Experimental ' yl]); ylabel(['Predicted ' yl]);
        if l == set_no
            legend('Experimental Data', ['Trained Data by GPtree ' num2str(l)])
        else
            legend('Experimental Data', ['Predicted Data by GPtree ' num2str(set_no) ' on dataset ' num2str(l)])
        end
        plot(DataSet(traindata1,f_index), traindata2, '.r');
        if saveon
            saveas(fig(2), [savedir '\' yl '_fit' num2str(set_no) 'on' num2str(l)], 'jpg')
            saveas(fig(1), [savedir '\' yl '_GPtree' num2str(set_no) 'on' num2str(l)], 'jpg')
            saveas(fig(2), [savedir '\' yl '_fit' num2str(set_no) 'on' num2str(l)], 'fig')
            saveas(fig(1), [savedir '\' yl '_GPtree' num2str(set_no) 'on' num2str(l)], 'fig')
        end
        pause(1)
    end
end
for i = 1:length(plst)
    rmpath([pwd plst{i}]);
end
end