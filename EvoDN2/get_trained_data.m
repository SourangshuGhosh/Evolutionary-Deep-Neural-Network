function mico_table = get_trained_data(parameters, set_no, yl,savedir)

global figure_han saveon ploton
setno = 1;
select = 1;
Xmin = parameters.Xmin;
Xmax = parameters.Xmax;
nooutnodes = parameters.nooutnodes;
f_index = parameters.out_index;
DataSet = parameters.DataSet;
ymin = parameters.Data_min(1,f_index);
ymax = parameters.Data_max(1,f_index);
Pop_str = parameters.Pop_str;
net = (parameters.dataset.pareto.P(select));
endnet = (parameters.dataset.pareto.P(select).endnet);
num_subnets = length(Pop_str);
in_end = [];
%select = parameters.dataset.pareto.select;
fig = figure_han;

mico_table = zeros(1,parameters.no_run);
for l = 1:parameters.no_run
    
    in = parameters.dataset(setno).in;
    out = parameters.dataset(setno).out;
    outmax = max(out);
    outmin = min(out);
    no_points = length(out);
    for subnet = 1:num_subnets
        in = parameters.dataset(setno).in;
        in = in(:,Pop_str{subnet}{1});
        no_layer = length(Pop_str{subnet}{2});
        snet = net.subnet{subnet};
        for layer = 1:no_layer - 1  %-1?
            in = in*snet{layer}(2:end,:);
            bias = snet{layer}(1,:);
            d=size(in);
            for i=1:d(1)
                for j=1:d(2)
                    in(i,j) = in(i,j)+bias(j);
                end
            end
            in = 1./(1+exp(-in));
        end
        in_end = [in_end in];
    end
    hold on
    modelout = in_end*endnet;
    % rescale the outputs
%     for i = 1:nooutnodes
%         modelout(:,i) = ymin(i)+(ymax(i)-ymin(i))*(modelout(:,i)-Xmin)/(Xmax-Xmin);
%         out(:,i) = ymin(i)+(ymax(i)-ymin(i))*(out(:,i)-Xmin)/(Xmax-Xmin);
%     end
    noexp = length(out(:,1));
    mico_table(1,l) = sqrt(sum(((modelout-out)).^2)/noexp);
    
    traindata1 = parameters.dataset(l).data_index;%DataSet(parameters.dataset(l).data_index, f1_index);
    traindata2 = modelout;
    if l == set_no && l == parameters.no_run, ploton = 1; end
    if ploton
        set(0,'CurrentFigure',fig(1)); clf
        plot(traindata1, DataSet(traindata1,f_index), '-rd'); hold on
        if l == set_no
            plot(traindata1, traindata2, '--ks')
            legend('Experimental Data', ['Trained Data by NNet ' num2str(l)])
        else
            plot(traindata1, traindata2, '--cs')
            legend('Experimental Data', ['Predicted Data by NNet ' num2str(set_no) ' on dataset ' num2str(l)])
        end
        xlabel('data index')
        ylabel(yl)
        pause(1)
        set(0,'CurrentFigure',fig(2)); hold off
        plot(DataSet(traindata1,f_index), DataSet(traindata1,f_index), '-kd', 'LineWidth', 2); hold on
        pol = polyfit(DataSet(traindata1,f_index), traindata2, 1);
        text(0.5*(0.5*(ymin+ymax)+ymax), 0.5*(0.5*(ymin+ymax)+ymin), ['slope of fitted line: ' num2str(pol(1))], 'Fontsize', 12)
        pol = pol(1)*DataSet(traindata1,f_index) + pol(2);
        plot(DataSet(traindata1,f_index), pol, '.-b', 'LineWidth', 2);
        xlabel(['Experimental ' yl]); ylabel(['Predicted ' yl]);
        if l == set_no
            legend('Experimental Data', ['Trained Data by NNet ' num2str(l)])
        else
            legend('Experimental Data', ['Predicted Data by NNet ' num2str(set_no) ' on dataset ' num2str(l)])
        end
        plot(DataSet(traindata1,f_index), traindata2, '.r');
        if saveon
            saveas(fig(2), [savedir '/' yl '_fit' num2str(set_no) 'on' num2str(l)], 'jpg')
            saveas(fig(1), [savedir '/' yl '_NNet' num2str(set_no) 'on' num2str(l)], 'jpg')
            saveas(fig(2), [savedir '/' yl '_fit' num2str(set_no) 'on' num2str(l)], 'fig')
            saveas(fig(1), [savedir '/' yl '_NNet' num2str(set_no) 'on' num2str(l)], 'fig')
        end
        pause(1)
    end
end

end