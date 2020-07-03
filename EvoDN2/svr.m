function response = svr(Setslog, figl, trend)

%Setslog = importdata(Setslog);
in_index = Setslog.in_index;
trend = importdata(trend);
response = [];
Xmin = Setslog.Xmin;
Xmax = Setslog.Xmax;
avg = ones(1,length(Setslog.in_index))*(Xmin+Xmax)/2;
for k = 1:Setslog.no_run
    setno = 1;
    select = 1;
    noinnodes = Setslog.noinnodes;
    %nonodes = Setslog.nonodes;
    nooutnodes = Setslog.nooutnodes;
    Pop_str = Setslog.Pop_str;
    net = (Setslog.dataset.pareto.P(select));
    endnet = (Setslog.dataset.pareto.P(select).endnet);
    num_subnets = length(Pop_str);
    in_end = [];
    %select = Setslog.dataset.pareto.select;
    for l = 1:noinnodes
         Prey = ones(length(trend),1)*avg;
         Prey(:,l)= trend;
         noexp = length(Prey(:,1));
         in_end = [];
         for subnet = 1:num_subnets
             in = Setslog.dataset(setno).in;
             in = Prey(:,Pop_str{subnet}{1});
             no_layer = length(Pop_str{subnet}{2});
             snet = net.subnet{subnet};
             for layer = 1:no_layer-1  %-1?
                 
                 s = in*snet{layer}(2:end,:);
                 bias = snet{layer}(1,:);
                 d=size(s);
                 for i=1:d(1)
                     for j=1:d(2)
                         s(i,j) = s(i,j)+bias(j);
                     end
                 end
                 in = 1./(1+exp(-s));
             end
             in_end = [in_end in];
             
         end
         hold on
         out = in_end*endnet;
         for i=1:nooutnodes
            if min(out) == max(out)
                out = 0.5*ones(size(out));
            else
                out = (out-min(out))/(max(out)-min(out));
            end
        end
        
        if k == Setslog.no_run
            hold off
            plot(1:length(Prey(:,1)), Prey(:,l), 'black','linewidth',4); hold on
            plot(1:length(Prey(:,1)), out, 'blue','linewidth',4);
            legend(gca,'  Input ', '  Output ');
            %legend(gca,'Nucleation density (input)', 'Grain size (output)');
            set(gca,'FontSize',20);
            scrsz = get(0, 'Screensize');
            set(gcf, 'OuterPosition', scrsz, 'Color', 'w')
        end
        p=diff(out);
        q=diff(trend);
        r=p.*q;
        maxr=max(r);
        minr=min(r);
        if(maxr<=0 && minr<=0)
            response(k,l)=-1;
            s = '-ve';
            %strep = 'Dirctly inverse';
        end
        if((maxr>=0 && minr>=0))
            response(k,l)=1;
            s = '+ve';
            %strep = 'Directly proportional';
        end
        if (maxr==0 && minr==0)
            response(k,l)=0;
            s = 'nil';
            %strep = 'No response';
        end
        if(minr<0 && maxr>0)
            response(k,l)=2;
            s = 'mix';
            %strep = 'Mixed response';
        end
        %title({['Testing parameter-' num2str(abc)]; strep});
        
        if k == Setslog.no_run
            text(270, 0.5, s, 'fontsize', 24)
            s = [ figl '-' num2str(k) '-' num2str(in_index(l))];
            %         title(gca, s,'fontsize',18);
            saveas(gcf, s, 'jpg');
            saveas(gcf, s, 'fig');
            pause(0.5)
        end
    end
    
end

disp(response)
save svr_res response

end
         
         