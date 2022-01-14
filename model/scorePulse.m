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

function [J, varargout] = scorePulse(x, SeasonLength, Seasons, m, p, IC, S, Method, sendToWorkspace)
%% Set defaults for optional arguments

    % rows is the number of potential pulses to score
    % cells is the number of plots
    % N is the total number of time intervals
    [rows, dummy] = size(x);
    [dummy, cells]= size(IC);
    N = Seasons*SeasonLength;
    have_MEX = exist('growAB_RK4');
    STEPS = 5;
    
    % set default values for optional arguments
    if(nargin<9)
        sendToWorkspace = 0;
        if(nargin<8)
            if(cells==1)
                Method =0;
            else
                Method =1;
            end
            
            if(nargin<7)
                S = 1;
            end
        end
    end
 
%% Set up Problem
    
    %Initialize Above and Below Ground 
    A = zeros(rows, cells); 
    if(length(IC(:,1)) == 1)
        B = repmat(IC, rows, 1);
    else
        B = IC;
    end
    
    if(length(m) > 1)
        m = repmat(m, 1, cells);
        p = repmat(p, 1, cells);
    end
    
    if(sendToWorkspace || Method > 1)
        As = zeros(rows, cells, N);
    end
    if(sendToWorkspace || Method ==1)
        Bs = zeros(rows, cells, N);
    end
    %Build Goat Deployment Matrix to Index into A
    % G is a logical matrix of size (rows, cells, N)
    % G1 is a temporary logical matrix of size (rows, cells*N) into which
    %   x can index directly.
    % c is a vector of integers 1:N, each repeating cells times.
    %  e.g. with cells=3 : [1 1 1 2 2 2 3 3 3 4 4 4 ... N N N]
    %  which indexes into G1
    G = false(rows, N*cells);
    for i=1:rows
        G(i, x(i,:)) = true;
    end
    G = reshape(G, [rows, cells, N]);
    
    %Find the Runge Kutta Time Step (h), do some precalculations.
    h = 1/(STEPS*SeasonLength);
    h2 = h/2;
    h6 = h/6;
    t = 0;
    S = (eye(size(S))+S*h);

%% Solve ODE System
    %Loop through each time interval
    for j = 1:N
        % We use an RK time step 1/5 of each time interval for precision
        for i = 1:STEPS
           
            %Runge Kutta-4 Method
            if(have_MEX)
                [A B] = growAB_RK4(A,B,h,m,p);
                t = t + h;
            else
                [k1A k1B] = growAB(A, B, t);
                [k2A k2B] = growAB(A+h2.*k1A, B+h2.*k1B, t);
                [k3A k3B] = growAB(A+h2.*k2A, B+h2.*k2B, t);
                [k4A k4B] = growAB(A+h.*k3A, B+h.*k3B, t);
                A = A + h6 * (k1A + 2*k2A + 2*k3A + k4A);  
                B = B + h6 * (k1B + 2*k2B + 2*k3B + k4B);
                t = t + h;
            end
            
            %Spread Above-Ground. This belongs propery in the growAB
            %  function but is placed here for speed.  This is equiv. to
            %  using Euler's Method (or RK1) for the spread portion.
            A = A*S;
        end
        
        %Store the results of this interval
        if(sendToWorkspace || Method > 1)
            As(:,:,j)=A;
        end
        if(sendToWorkspace || Method ==1)
            Bs(:,:,j)=B;
        end
        
        % If it's the end of a season, set Above-Ground to 0
        %   note that t is incremented by h in the nested loop below
        %   and has units of seasons
        if(mod(j, SeasonLength)==0)
            A = zeros(rows, cells);
        end
        
        % Set Above-Ground to 0 if there's a goat deployment this interval
        A(G(:,:,j)) = 0;
  
    end
    
%% Find Value of the Objective Function     
    if(Method == 0)
        % Below-Ground of cell 1 at the end
        J = B(:,1);
    elseif(Method ==1)
        % Summed Below-Ground of all cells at the end
        J = sum(B,2);
    elseif(Method ==2)
        % Summed Above Ground of All cells through all time
        J = sum(sum(As,3),2);
    elseif(Method ==3)
        % Summed Above Ground of Central cell through all time
        J = sum(As(:,1,:),3);
    end
    
%% Place some values in the Workspace if Requested
    if(sendToWorkspace)
        varargout(1) = {permute(As, [3 2 1])};
        varargout(2) = {permute(Bs, [3 2 1])};
        varargout(3) = {permute(G, [3 2 1])};
    end
    
%% Define the Growth Function, without Spread
    function [dA dB] = growAB(A, B, t)
        dA = m.*B.*(1-A);
        dB = p.*A - dA;
    end
end