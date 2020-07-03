function response = svr(filename, trend)

Setslog = importdata(filename);
in_index = Setslog.in_index;
trend = importdata(trend);
response = [];
Xmin = Setslog.Xmin;
Xmax = Setslog.Xmax;
avg = ones(1,length(Setslog.in_index))*(Xmin+Xmax)/2;
for k = 1:Setslog.no_run
%     in = Setslog.dataset(k).in;
    
    noinnodes = Setslog.noinnodes;
    nonodes = Setslog.nonodes;
    nooutnodes = Setslog.nooutnodes;
    i = Setslog.dataset(k).pareto.select;
    Net = Setslog.dataset(k).pareto.P(i);
    w = Net.w;
    W = Net.W;
    for l = 1:noinnodes
        Prey = ones(length(trend),1)*avg;
        Prey(:,l)= trend;
        
        noexp = length(Prey(:,1));
        s = zeros(noexp,1); z = [];
        for i = 1:nonodes
            s = w(i,1)*ones(noexp,1)+((w(i,2:noinnodes+1)*Prey')');
            z(:,i) = 1./(1+exp(-s));
        end
        
        A = [ones(noexp,1) z];
        
        out = A*W;
        
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
            legend(gca,'  Input ', '  Output ','FontSize',20);
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

disp(response)
save svr_res response

end