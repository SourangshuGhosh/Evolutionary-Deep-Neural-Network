%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  The source code of the reference vector guided evolutionary algorithm (RVEA)
%%
%%  See the details of RVEA in the following paper:
%%
%%  R. Cheng, Y. Jin, M. Olhofer and B. Sendhoff, 
%%  A Reference Vector Guided Evolutionary Algorithm for Many-objective Optimization,
%%  IEEE Transactions on Evolutionary Computation, 2016
%%
%%  The source code of RVEA is implemented by Ran Cheng 
%%
%%  If you have any questions about the code, please contact: 
%%  
%%  Ran Cheng at ranchengcn@gmail.com
%%  Prof. Yaochu Jin at yaochu.jin@surrey.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Opt(Problem,TrainType,parameters,savedir)
oldpath = path;
clc;format short;
newpath = pwd;
addpath([newpath '\Public'], [newpath '\RVEA']);
warning off
Objectives = length(parameters.out_index);
Algorithm = {'RVEA'};
MAIN(TrainType,Objectives,parameters,savedir)
path(oldpath);
end