function con_val = Contraints(Algorithm,Var,varargin)
%contra - Constraint handling in BioGP, EvoNN and cRVEA
% Version 0.5
% Syntax: con_val = contra(input)
%
% Long description
	switch Algorithm
	case 'BioGP'
		F1 = varargin{1};
		F2 = varargin{2};
		con_val = [F1 F2];
		con_val = contra(Var,F1,F2);
	case 'EvoNN'
		F1 = varargin{1};
		F2 = varargin{2};
		con_val = [F1 F2];
		con_val = contra(Var,F1,F2);
	case 'cRVEA'
		Obj = varargin{1};
		eqCon = varargin{2};
		ieqCon = varargin{3};
        obj = varargin{4};
		con_val = zeros(size(Obj,1));
		con_val = contra_cRVEA(Var,Obj,eqCon,ieqCon,obj);
	end
end

function Fout = contra(Prey, F1, F2)  % constraint function
	global obj LB_F UB_F
	%penalty = 1e1;
	F = [F1 F2];
	flag = 0;
	for i = 1:length(F1)
		for j = 1:length(F(1,:))
			if F(i,j) < LB_F(j) || F(i,j) > UB_F(j)
				flag = 1;
			end
		end
		if flag
			F1(i) = F1(i)*(1 +obj(1)*0.5);       
			F2(i) = F2(i)* (1 + obj(2)*0.5);
		end
		flag = 0;
        Fout = [F1 F2];
	end
end

function con_val = contra_cRVEA(Var, Obj, eqCon, ieqCon,obj)
	numVar = size(Var,2);
	numObj = size(Obj,2);
	numEC = length(eqCon);
	numIC = length(ieqCon);
    eqSum = zeros((size(Obj,1)),1);
    ieqSum = eqSum;
	if numEC == 1
		if length(eqCon{1}) == 0
			numEC = 0;
		end
	end
		if numIC == 1
		if length(ieqCon{1}) == 0
			numIC = 0;
		end
	end
	for i = 1:numVar
		eval(['var' num2str(i) ' = Var(:,' num2str(i) ');']);
	end
	for i = 1:numObj
		eval(['obj' num2str(i) ' = obj(' num2str(i) ')*Obj(:,' num2str(i) ');']);
    end
	for i = 1:numEC
		eval(['eqVal' num2str(i) ' = ' eqCon{i} ';']);
		eval(['eqVal' num2str(i) '(eqVal' num2str(i) '<0) = -eqVal' num2str(i) '(eqVal' num2str(i) '<0);']);
		eval(['eqSum = eqSum + eqVal' num2str(i) ';']);
	end
	for i = 1:numIC
        %disp(['ieqVal' num2str(i) ' = ' ieqCon{i} ';'])
		eval(['ieqVal' num2str(i) ' = ' ieqCon{i} ';']);
		eval(['ieqVal' num2str(i) '(ieqVal' num2str(i) '>0) = 0;']);
		eval(['ieqVal' num2str(i) '(ieqVal' num2str(i) '<0) = -ieqVal' num2str(i) '(ieqVal' num2str(i) '<0);']);
        eval(['ieqSum = ieqSum + ieqVal' num2str(i) ';']);
	end
	con_val = eqSum + ieqSum;
end