function response = svr_BioGP(filename)
cd ..\BioGP
plst = {'\adf'};
for i = 1:length(plst)
    path([pwd plst{i}], path);
end
cd ..\Tools
Setslog = importdata([filename '.mat']);
T1_set = Setslog.T.set;
in_index = Setslog.in_index;
trend = importdata('trend.mat');
response = zeros(Setslog.no_run, length(in_index));
in = Setslog.dataset(Setslog.no_run).in;
mx = max(in);
mn = min(in); avg = (mx + mn)/2;
for k = 1:Setslog.no_run
    
    p = Setslog.dataset(k).pareto.P(1);
    
    for l = 1:length(in(1,:))
        Prey = ones(length(trend),1)*avg;
        Prey(:,l)= mn(l) + (mx(l)-mn(l))*trend;
        
        for i = 1:length(in_index)
            eval([T1_set{i} ' = Prey(:,i);']);
        end
        noexp = length(Prey(:,l));
        feval = ones(noexp, p.no_roots+1); w = ones(p.no_roots+1,1); w(1) = p.bias;
        for i = 1:p.no_roots
            tree = eval(cat(2,p.root(i).tree{:}));
            feval(:,i+1) = tree.*ones(noexp,1);
            w(i+1) = p.root(i).w;
        end
        out = feval*w;
        
        if min(out) == max(out)
            out = 0.5*ones(size(out));
        else
            out = (out-min(out))/(max(out)-min(out));
        end
        
        if k == Setslog.no_run
            hold off
            plot(1:length(Prey(:,1)), trend, 'g','linewidth',2); hold on
            plot(1:length(Prey(:,1)), out, '--r','linewidth',2);
            legend(gca,'Input Signal', 'Response Signal');
            scrsz = get(0, 'Screensize');
            set(gcf, 'OuterPosition', scrsz, 'Color', 'w')
        end
        r=diff(out);
        q=diff(trend);
        r=r.*q;
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
            text(length(trend)+20, 0.5, s, 'fontsize', 30)
            foldername = [filename '_SVR'];
            mkdir(foldername);
            s = [foldername '\run_no' num2str(k) '-var' num2str(in_index(l))];
            %         title(gca, s,'fontsize',18);
            saveas(gcf, s, 'jpg');
            saveas(gcf, s, 'fig');
            pause(0.5)
        end
    end
    
end
cd ..\BioGP
for i = 1:length(plst)
    rmpath([pwd plst{i}]);
end
cd ..\Tools
disp(response)
save svr_res2 response
end