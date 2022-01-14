% [x,fval,exitflag,output,population,score] = ChooseDaysGA(pulses, cells, SeasonLength, Seasons, m, p, IC, S, Method,ShowOutput)
% pulses: the length of x; the number of total pulses over all cells and
%        all seasons
% cells: the number of cells for multi-cell problems
% SeasonLength: the number of intervals within a season
% Seasons: the total number of seasons
% m, p: model parameters
% IC: vector of initial conditions for B0 for each cell
% S: a spread matrix for multi-cell problems
% Method: the scoring method.  0 for B(T) in cell 1, 1 for sum(B(T)) over
%    all cells, and 2 for sum(A) over all time intervals and cells.
% ShowOutput: 0 to repress output, 1 to watch progress.

function [x,fval,exitflag,output,population,score] = chooseDaysGA(pulses, cells, SeasonLength, Seasons, m, p, IC, S, Method,ShowOutput)

    % Start with the default options
    options = gaoptimset;

    % Show or Repress Output.  Always repress when running on the cluster
    if(nargin<10)
        ShowOutput =1;
    end
    if(ShowOutput ==1)
        options = gaoptimset(options,'Display', 'iter');
    else
        options = gaoptimset(options,'Display', 'none');
    end
    
    % Set the upper bound of values on X to pass to scorePulses()
    ubound = Seasons*SeasonLength*cells;
    
    % Modify options setting
    options = gaoptimset(options,'PopInitRange', [1;ubound]);
    options = gaoptimset(options,'PopulationType', 'custom');
    options = gaoptimset(options,'CreationFcn', @init_pop, 'CrossoverFcn', @crossoverSingle,  'MutationFcn', {@mutationInteger, 0.2});
    options = gaoptimset(options,'PopulationSize', ubound*10);
    options = gaoptimset(options,'UseParallel', 'never');
    options = gaoptimset(options,'Vectorize', 'on');
    options = gaoptimset(options,'EliteCount', 10);
    options = gaoptimset(options,'CrossoverFraction', 0.7);
    options = gaoptimset(options,'StallGenLimit', 5, 'Generations', 100, 'TimeLimit', 900, 'TolFun', 1e-7);
    %options = gaoptimset(options,'SelectionFcn', @selectionChoose);
    %options = gaoptimset(options,'FitnessScalingFcn', @fitnessrank);
    options = gaoptimset(options,'OutputFcns', { [] });

    [x,fval,exitflag,output,population,score] = ga({@scorePulse, SeasonLength, Seasons, m, p, IC, S, Method},pulses,[],[],[],[],1,ubound,[],options);


% Population Creation Fcn
function Population = init_pop(GenomeLength,FitnessFcn,options)
    lower = options.PopInitRange(1,end);
    upper = options.PopInitRange(2,end);
    Population = randi([lower upper], options.PopulationSize, GenomeLength);
    
  
% Mutation Function
function kids = mutationInteger(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation, MutationRate)
    lower = options.PopInitRange(1,end);
    upper = options.PopInitRange(2,end);
   
    kids = thisPopulation(parents,:);
   
    mutations = (rand(length(parents), GenomeLength) < MutationRate);
    newVals = randi([lower upper], sum(sum(mutations)),1);
    kids(mutations) = newVals;
    
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