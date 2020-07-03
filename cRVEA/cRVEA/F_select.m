function [Selection] = F_select(FunctionValue, V, theta0, refV, FinalSelection)

CV = FunctionValue(:,end);
FunctionValue = FunctionValue(:,1:end - 1);
[N M] = size(FunctionValue);
VN = size(V, 1);
%Function Value Translation
Zmin = min(FunctionValue,[],1);
FunctionValue = (FunctionValue - repmat(Zmin, [size(FunctionValue,1) 1]));

%Solutions associattion to reference vectors
clear class;
uFunctionValue = FunctionValue./repmat(sqrt(sum(FunctionValue.^2,2)), [1 M]);
cosine = uFunctionValue*V'; %calculate the cosine values between each solution and each vector
acosine = acos(cosine);
[maxc maxcidx] = max(cosine, [], 2);
class = struct('c', cell(1,VN)); %classification
for i = 1:N
    class(maxcidx(i)).c = [class(maxcidx(i)).c, i];
end;

Selection = [];
for k = 1:VN
    if(~isempty(class(k).c))
        sub = class(k).c;
        subFunctionValue = FunctionValue(sub,:);
        
        %APD calculation
        subacosine = acosine(sub, k);
        subacosine = subacosine/refV(k);
        D1 = sqrt(sum(subFunctionValue.^2,2));
        D = D1.*(1 + (theta0)*(subacosine));
        
        %constraints violence based selection
        subCV = CV(sub, :);
        fsSet = find(subCV == 0);
        if(length(fsSet)) > 0
            %if there exists any feasible solution, select the one with the
            %smalles APD
            D = D(fsSet);
            [mind mindidx] = min(D);
            sub = sub(fsSet);
            Selection = [Selection; sub(mindidx)];
        else
            %if all solutions are infeasible, and it is not the final selection, select the one with the smallest
            %violence
            if(~FinalSelection)
                [mind mindidx] = min(subCV);
                Selection = [Selection; sub(mindidx)];
            end;
        end;
        
    end;
    
end;


end

