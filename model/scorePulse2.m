%[J, varargout] = scorePulse(x, SeasonLength, Seasons, m, p, IC, S, Method, sendToWorkspace)
% x is a vector with values specifying when to deploy goats.
%   values equal Time*(NumberOfCells-1)+CellNumber
% SeasonLength is the number of time intervals in a single season
% Seasons is the total number of seasons
% m, p are model parameters
% IC is a vector of initial conditions on B per cell
% S is an optional Spread Matrix for multi-plot problems
% Method: the scoring method.  0 for B(T) in cell 1, 1 for sum(B(T)) over
%    all cells, and 2 for sum(A) over all time intervals and cells.
% sendToWorkspace is true to return G, A, and B along with the score, J

function [J, varargout] = scorePulse(x, SeasonLength, Seasons, m, p, IC, Method, sendToWorkspace)
%% Set up Problem
    %Initialize Above and Below Ground 
    N = Seasons*SeasonLength;
    A = 0; 
    B = IC;
    As = zeros(N,1);
    Bs = zeros(N,1);
   
    G = zeros(N,1);
    G(x) = 1;
    
    %Find the Runge Kutta Time Step (h), do some precalculations.
    h = 1/(5*SeasonLength);
    h2 = h/2;
    h6 = h/6;
    epsilon = 1/(10*N);
    t = 0;

%% Solve ODE System
    %Loop through each time interval
    for j = 1:N
        if(t-floor(t) < epsilon)
            A = 0;
        end
        if G(j),
            A = 0;
        end
        
        for i = 1:5
            %Runge Kutta-4 Method
            [k1A k1B] = growAB(A, B, t);
            [k2A k2B] = growAB(A+h2.*k1A, B+h2.*k1B, t);
            [k3A k3B] = growAB(A+h2.*k2A, B+h2.*k2B, t);
            [k4A k4B] = growAB(A+h.*k3A, B+h.*k3B, t);
            A = A + h6 * (k1A + 2*k2A + 2*k3A + k4A);  
            B = B + h6 * (k1B + 2*k2B + 2*k3B + k4B);
            t = t + h;
        end
        
        As(j)=A;
        Bs(j)=B;
    end
    
%% Find Value of the Objective Function     
    J = [B sum(A)];
    
%% Place some values in the Workspace if Requested
    if(sendToWorkspace)
        varargout(1) = {permute(As, [3 2 1])};
        varargout(2) = {permute(Bs, [3 2 1])};
        varargout(3) = {permute(G, [3 2 1])};
    end
    
    %% Define the Growth Function, without Spread
    function [dA dB] = growAB(A, B, t)
        dA = m*B.*(1-A);
        dB = p*A - dA;
    end
    
end