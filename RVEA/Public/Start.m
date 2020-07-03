function Start(Algorithm,Problem,Objectives,Run)
    if nargin < 4;Run = 0;end
    if ~iscell(Algorithm);Algorithm = {Algorithm};end
    if ~iscell(Problem);Problem = {Problem};end
    for R = Run
        for A = Algorithm
            a = cell2mat(A);
            if exist(a,'dir') == 7
                addpath(a);
            else
                error(['Algorithm ',a,' does not exist.']);
            end
            for P = Problem
                for M = Objectives
                    MAIN(cell2mat(P),M,R);
                end
            end
        end
    end
    %sp=actxserver('SAPI.SpVoice');
    %sp.Speak('job finished');
end

