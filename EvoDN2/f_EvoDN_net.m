function f_EvoDN_net(parameters)
    setno = 1;
    select = 1;
	in = parameters.dataset(setno).in;
	out = parameters.dataset(setno).out;
    Pop_str = parameters.Pop_str;
    net = (parameters.dataset.pareto.P(select));
    endnet = (parameters.dataset.pareto.P(select).endnet);
    num_subnets = length(Pop_str);
    in_end = [];
    select = parameters.dataset.pareto.select;
	
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
            in = in+bias;
            in = 1./(1+exp(-in));
        end
        in_end = [in_end in];
    end
    hold on
 	modelout = in_end*endnet;
    plot(out,modelout,'kd','LineWidth', 2)
    pol = polyfit(out,modelout,1);
    ylim=get(gca,'ylim');
    xlim=get(gca,'xlim');
    text(xlim(1)+0.1*(xlim(2)-xlim(1)),ylim(2),['\fontsize{12} \color{red} Slope =' num2str(pol(1))])
    pol = pol(1)*out + pol(2);
    plot(out,pol, '.-b', 'LineWidth', 2)
    %plot([1:no_points],out,'--o',[1:no_points],modelout,'-or')
    hold off

end