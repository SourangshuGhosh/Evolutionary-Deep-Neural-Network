function P_output (Population,time,Algorithm,Problem,M, obj_val,savedir,obj,eqCon,ieqCon)

FunctionValue = P_objective('value',obj_val,M,Population,obj,eqCon,ieqCon);
if(strcmp(Algorithm, 'cRVEA'))
    FunctionValue = FunctionValue(:,1:end - 1);
end;

NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue = FunctionValue(NonDominated,:);

if(M == 2)
    plot(obj(1)*FunctionValue(:,1), obj(2)*FunctionValue(:,2), 'ro', 'MarkerFace', 'r');
    %xlim([0 max(FunctionValue(:,1)) + 0.1]);
    %ylim([0 max(FunctionValue(:,2)) + 0.1]);
    xlabel('F1');ylabel('F2');
    hold off;
else
    plot3(obj(1)*FunctionValue(:,1), obj(2)*FunctionValue(:,2), obj(3)*FunctionValue(:,3), 'ro','MarkerFace', 'r');
%     xlim([0 max(FunctionValue(:,1)) + 0.1]);
%     ylim([0 max(FunctionValue(:,2)) + 0.1]);
%     zlim([0 max(FunctionValue(:,3)) + 0.1]);
    xlabel('F1', 'FontSize', 14);ylabel('F2', 'FontSize', 14);zlabel('F2', 'FontSize', 14);
    view(135, 30);
    hold off;
end;
eval(['save ' savedir '/cRVEAopt.mat Population FunctionValue time'])
saveas(gcf, [savedir '/cRVEA_Pareto'], 'jpg')
saveas(gcf, [savedir '/cRVEA_Pareto'], 'fig')
end


