function [x,fval,exitflag,output,population,score] = goatdaysoptim(nvars)
% This is an auto generated MATLAB file from Optimization Tool.


lb = 1;
ub = 98;
PopInitRange_Data = [1;98];


% Start with the default options
options = gaoptimset;
% Modify options setting
options = gaoptimset(options,'PopInitRange', PopInitRange_Data);
options = gaoptimset(options,'FitnessLimit', 5);
options = gaoptimset(options,'CrossoverFcn', {  @crossoverintermediate [] });
options = gaoptimset(options,'Display', 'final');
options = gaoptimset(options,'OutputFcns', { [] });
[x,fval,exitflag,output,population,score] = ...
gamultiobj(@chooseGoatDays,nvars,[],[],[],[],lb,ub,options);
