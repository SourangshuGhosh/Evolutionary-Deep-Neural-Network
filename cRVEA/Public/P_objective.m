function [Output, Boundary, Coding] = P_objective(Operation,wlog,M,Input,varargin)

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
            FunctionValue = [];
			obj = varargin{1};
			eqCon = varargin{2};
			ieqCon = varargin{3};
			for i = 1:M 
				FunctionValue = [FunctionValue obj(i)*evaluate_obj(Input, wlog(i))];
			end
            FunctionValue = [FunctionValue Contraints('cRVEA',Input,FunctionValue,eqCon,ieqCon,obj)];
			Output = FunctionValue;
    end
end