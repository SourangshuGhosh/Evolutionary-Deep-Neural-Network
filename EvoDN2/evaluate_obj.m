function out=evaluate_obj(data,parameters)
	%nonodes = parameters.nonodes;
    noinnodes = parameters.noinnodes;
    nooutnodes = parameters.nooutnodes;
    xmin = parameters.xmin;
    xmax = parameters.xmax;
    ymin = parameters.ymin;
    ymax = parameters.ymax;
    Pop_str = parameters.Pop_str;
    net = parameters.net;
    endnet = parameters.endnet;
    num_subnets = length(Pop_str);
    in_end = [];
    x= data(:,1:noinnodes);

    %scale the inputs
    Xmin = parameters.Xmin;
    Xmax = parameters.Xmax;
    xsc=[];
    % xsc=parameters.DataSet_sc;
    for i=1:noinnodes
        xsc=[xsc Xmin+ (x(:,i)-xmin(i))/(xmax(i)-xmin(i))*(Xmax-Xmin)];
    end
    in = xsc;
    
    noexp = length(in(:,1));
	
	for subnet = 1:num_subnets
        s = parameters.in;
        s = in(:,Pop_str{subnet}{1});
        no_layer = length(Pop_str{subnet}{2});
        snet = net.subnet{subnet};
        for layer = 1:no_layer - 1  %-1?
            s = s*snet{layer}(2:end,:);
            bias = snet{layer}(1,:);
            d=size(s);
                 for i=1:d(1)
                     for j=1:d(2)
                         s(i,j) = s(i,j)+bias(j);
                     end
                 end
            s = 1./(1+exp(-s));
        end
        in_end = [in_end s];
    end
    hold on
 	out = in_end*endnet;
end