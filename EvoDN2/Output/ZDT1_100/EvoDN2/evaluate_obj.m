function out=evaluate_obj(data,parameters)
	 Pop_str = parameters.Pop_str;
    net = (parameters.dataset.pareto.P(select));
    endnet = (parameters.dataset.pareto.P(select).endnet);
    num_subnets = length(Pop_str);
    in_end = [];
	in = data;
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
 	out = in_end*endnet;
end