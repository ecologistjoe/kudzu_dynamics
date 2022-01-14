function [J varargout] = chooseGoatDays2d(x, z, c, IC, sendToWorkspace)
%% Determine the Form of inputs and Create the Control Regime
    if(nargin<5)
        sendToWorkspace = 0;
    end

    [rows cols] = size(x);
    

    

 %% Set up Problem
    
    A = ones(rows,1) * IC(1);
    B = ones(rows,1) * IC(2);
    N = 200;
    h = 1/(5*N);
    t = h;
    h2 = h/2;
    h6 = h/6;
    Bs = zeros(rows, N);
    As = zeros(rows, N);
 
%% Solve ODE System
    for j = 1:N
        Bs(:,j)=B;
        As(:,j)=A;
        for i = 1:5
            [k1A k1B] = growAB(A, B, t);
            [k2A k2B] = growAB(A+h2.*k1A, B+h2.*k1B, t);
            [k3A k3B] = growAB(A+h2.*k2A, B+h2.*k2B, t);
            [k4A k4B] = growAB(A+h.*k3A, B+h.*k3B, t);
            A = A + h6 * (k1A + 2*k2A + 2*k3A + k4A);
            B = B + h6 * (k1B + 2*k2B + 2*k3B + k4B);
            t = t + h;
        end
    end
    
%% Find Value of the Objective Function     
    if(nargin<3)
        w = [3000 200 100 1];
    end
    
    R = 2*log(1+1./(5*B.^2)) ./ (m*(sqrt(5)-1)); % From Maple, see notes.
    o = zeros(rows,1);
    J = w(1)*(3-min(R,2) - sum(As<.1,2)/N) +...
        w(2)*sum(diff([o G o],1,2)~=0,2) + ...
        w(3)*sum(G>0,2) + ...
        w(4)*sum(max(G-100*As,0),2);
    
    J(J<0) = 1e10;
    J(isnan(J)) = 1e10;

%% Place some values in the Workspace if Requested
    if(sendToWorkspace)
        varargout(1) = {As'};
        varargout(2) = {Bs'};
        varargout(3) = {G'};
    end
    
%% Define the Growth Function
    function [dA dB] = growAB(A, B, t)
        dA = B.*(1-A) - G(:,j).*A .* (1.1)./(A+0.1);
        dB = z*A - B.*(1-A);
    end

end