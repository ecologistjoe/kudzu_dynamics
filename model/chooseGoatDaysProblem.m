function [x,fval,exitflag,output,population,score] = chooseGoatDaysProblem(nVars, weights, ICs)
% This is an auto generated MATLAB file from Optimization Tool.

% Start with the default options
options = gaoptimset;

if(nVars == 100)
    Problem = @chooseGoatDays;
elseif(nVars == 200)
    Problem = @chooseGoatDaysTwo;
elseif(nVars == 300)
    Problem = @chooseGoatDaysThree;
end

if(nVars >= 5 && nVars <=25)
    PopRange = [-1;1];
    options = gaoptimset(options,'PopulationType', 'doubleVector');
    options = gaoptimset(options,'MutationFcn', {@mutationDCT, .1});
    options = gaoptimset(options,'PopulationSize', 1000);
else
    if(nVars == 40)
        PopRange = [0;197];
    else
        PopRange = [0;10];
    end
    options = gaoptimset(options,'PopulationType', 'custom');
    options = gaoptimset(options,'CreationFcn', @init_pop, 'CrossoverFcn', @crossoverSingle,  'MutationFcn', {@mutationInteger, 0.2});
    options = gaoptimset(options,'PopulationSize', 1000+nVars*10);
end

% Modify options setting
options = gaoptimset(options,'UseParallel', 'Always');
options = gaoptimset(options,'Vectorize', 'on');
options = gaoptimset(options,'PopInitRange', [0;10]);
options = gaoptimset(options,'EliteCount', 10);
options = gaoptimset(options,'CrossoverFraction', 0.7);
options = gaoptimset(options,'StallGenLimit', 5, 'Generations', 100, 'TimeLimit', 900, 'TolFun', 1e-7);
options = gaoptimset(options,'SelectionFcn', @selectionChoose);
options = gaoptimset(options,'FitnessScalingFcn', @fitnessrank);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'OutputFcns', { [] });

[x,fval,exitflag,output,population,score] = ga({Problem,0, weights, ICs},nVars,[],[],[],[],0,10,[],options);


% Population Creation Fcn
function Population = init_pop(GenomeLength,FitnessFcn,options)
    range = options.PopInitRange;
    lower = range(1,:);
    upper = range(2,:);
    Population = zeros(options.PopulationSize, GenomeLength);
    
    Z = (rand(options.PopulationSize, GenomeLength) < .01);
    %Z = conv2(double(Z), [1 1 1], 'same');
    %Z = logical(Z);
    Population(Z) = 6;
    Population(1,:) = zeros(1,GenomeLength);
    Population(2,:) = ones(1,GenomeLength);
    
% End of creation function
%---------------------------------------------------


function expectation = fitnessrank(scores,nParents)
global bestScore bestMemberIndex;

scores = scores(:);
[J,i] = sort(scores);

expectation = zeros(size(scores));
expectation(i) = 1 ./ ((1:length(scores))  .^ 0.5);

expectation = nParents * expectation ./ sum(expectation);


% Mutation Function
function kids = mutationInteger(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation, MutationRate)

    lower = options.PopInitRange(1,1);
    upper = options.PopInitRange(2,1);
   
    range = round(upper-lower) - 0.5*(state.Generation/options.Generations);
    L = length(parents);
    kids = thisPopulation(parents,:);
   
    mutations = (rand(L, GenomeLength) < MutationRate);
    newVals = randi([0 round(range)], sum(sum(mutations)),1) - round(range/2);
    kids(mutations) = kids(mutations) + newVals;
    
    Z = (rand(L, GenomeLength) < MutationRate);
    kids(Z) = 0;
    
    %Bound
    kids(kids>upper) = upper;
    kids(kids<lower) = lower;

%---------------------------------------------------

% Mutation Function for the DCT selection
function kids = mutationDCT(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation, MutationRate) 
    kids = thisPopulation(parents,:);
    mutations = (rand(size(kids)) < MutationRate);
    newVals =  randn(sum(mutations(:)), 1);
    kids(mutations) = kids(mutations) + newVals;
    
    %Bound
    kids(kids>1) = 1;
    kids(kids<-1) = -1;
%---------------------------------------------------



%Selection Function
function parents = selectionChoose(expectations,nParents,options)
expectation = expectations(:,1);
wheel = cumsum(expectation) / nParents;

parents = ones(1,nParents);
% we will step through the wheel in even steps.
stepSize = 1/nParents;

% we will start at a random position less that one full step
position = rand * stepSize;

% a speed optimization. Position is monotonically rising.
lowest = 1; 

for i = 1:nParents % for each parent needed,
    for j = lowest:length(wheel) % find the wheel position
        if(position < wheel(j)) % that this step falls in.
            parents(i) = j;
            lowest = j;
            break;
        end
    end
    position = position + stepSize; % take the next step.
end

% End Selection Function


% Cross-Over Function
function kids  = crossoverSingle(parents,options,GenomeLength,FitnessFcn,unused,thisPopulation)
% How many children to produce?
nKids = length(parents)/2;
kids = zeros(nKids,GenomeLength);
x = randi(GenomeLength-1,1,nKids);    
    
for i=1:nKids
    % get parents
    parent1 = thisPopulation(parents(2*i-1),:);
    parent2 = thisPopulation(parents(2*i),:);
    % cut point is AFTER this index.
    kids(i,:) = [parent1(1:x(i)) parent2((x(i)+1):end)];
end