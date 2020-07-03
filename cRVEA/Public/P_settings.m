function [Generations,N,varargout] = P_settings(Algorithm,Problem,M,parameters)
%myFun - Description
%
% Syntax: [Generations,N,varargout] = myFun(input)
%
% Long description
	Generations = parameters.cRVEAopt.Generations;
	varargout = parameters.cRVEAopt.p1p2;
	N = parameters.cRVEAopt.N;
end