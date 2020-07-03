function [Output, Boundary, Coding] = P_objective(Operation,wlog,M,Input)

	Boundary = NaN; Coding = NaN;
	switch Operation
		case 'init'
			novar = length(wlog(1).in_index);
			LB = []; UB = [];
				for i = 1:M
				LB = [LB; wlog(i).xmin];
				UB = [UB; wlog(i).xmax];
			end
			LB = max(LB); UB = min(UB);
			Output = rand(Input,novar);
			for i = 1:novar
				Output(:,i) = Output(:,i)*(UB(i)-LB(i))+LB(i);
			end
			Coding = 'Real'; Boundary = [UB;LB];
		case 'value'
			Population = Input;
			for i = 1:M 
				FunctionValue(:,i) = evaluate_obj(Input, wlog(i));
			end
			Output = FunctionValue;
    end
end