function P_output (Population,time,Algorithm,Problem,M, obj_val,savedir)

FunctionValue = P_objective('value',obj_val,M,Population);
if(strcmp(Algorithm, 'cRVEA'))
    FunctionValue = FunctionValue(:,1:end - 1);
end;
%TrueValue = P_objective('true',Problem,M,1000);
%del = find(FunctionValue(:,1) < 0)
%FunctionValue(del,:) = [];
%del = find(FunctionValue(:,2) < 0)
%FunctionValue(del,:) = [];
NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue = FunctionValue(NonDominated,:);


if(M == 2)
	plot(obj(1)*FunctionValue(:,1), obj(2)*FunctionValue(:,2), '--ro');
    xlabel('f_1');ylabel('f_2');
else
    plot3(TrueValue(:,1), TrueValue(:,2), TrueValue(:,3), '.');
    hold on;
    plot3(obj(1)*FunctionValue(:,1), obj(2)*FunctionValue(:,2), obj(3)*FunctionValue(:,3), 'ro','MarkerFace', 'r');
%     xlim([0 max(TrueValue(:,1)) + 0.1]);
%     ylim([0 max(TrueValue(:,2)) + 0.1]);
%     zlim([0 max(TrueValue(:,3)) + 0.1]);
    xlabel('f_1', 'FontSize', 14);ylabel('f_2', 'FontSize', 14);zlabel('f_3', 'FontSize', 14);
    view(135, 30);
    hold off;
end;

eval(['save ' savedir '/RVEAopt.mat Population FunctionValue time'])
end


