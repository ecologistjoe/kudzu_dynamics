function [x,fval,exitflag,output,population,score] = goatoptim(nvars,lb,ub,PopInitRange_Data,PopulationSize_Data,StallGenLimit_Data)
% This is an auto generated MATLAB file from Optimization Tool.

% Start with the default options
options = gaoptimset;
% Modify options setting
options = gaoptimset(options,'PopInitRange', PopInitRange_Data);
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'StallGenLimit', StallGenLimit_Data);
options = gaoptimset(options,'FitnessScalingFcn', @fitscalingprop);
options = gaoptimset(options,'CrossoverFcn', @crossoverscattered);
options = gaoptimset(options,'MutationFcn', {  @mutationuniform 0.05 });
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'OutputFcns', { [] });
[x,fval,exitflag,output,population,score] = ...
ga(@chooseGoatDays,nvars,[],[],[],[],lb,ub,[],options);
x