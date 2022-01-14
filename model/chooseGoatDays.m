function [J varargout] = chooseGoatDays(x,sendToWorkspace, w, IC)
%% Determine the Form of inputs and Create the Control Regime
    if(nargin<2)
        sendToWorkspace = 0;
    end

    [rows cols] = size(x);
    
    if(cols == 200)
        G = x*10;
    elseif(cols ==100)
        G = x(:,floor(1:.5:100.5))*10;
    elseif(cols == 40)
        G = zeros(rows, 200);
        s = x(:,1:2:39)+1;
        g = x(:,2:2:40)/10;
        
        G(s) = g;
        G(s+1) = g;
        G(s+2) = g;
    elseif(cols == 50)
        %duplicate values to make twice as long
        G = x(:,floor(1:.25:50.75));
        cols = 50;
    elseif(cols>=5 && cols<=25)
        %X's are DCT values
        G = idct(x', 200)'+.2;
        G(G<0) = 0;
        G = round(G*10)*10;
    end
    
%Step function spacing
%     x = ceil(x);
%     if((sum(x<1) + sum(x>99)) > 0)
%         J = 1e10;
%         return;
%     end
%     
%     G = zeros(1,100);
%     G([x x+1 x+2]) = 1;
    

 %% Set up Problem
    
    A = ones(rows,1) * IC(1);
    B = ones(rows,1) * IC(2);
    m = 7;
    p = 7;
    N = 200;
    h = 1/(5*N);
    t = h;
    h2 = h/2;
    h6 = h/6;
    Bs = zeros(rows, N);
    As = zeros(rows, N);
 
%% Solve ODE System
    for j = 1:N
        if(j == 200)
            A = zeros(rows,1);
        end
        if(j == 400)
            A = zeros(rows,1);
        end
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
        dA = m*B.*(1-A) - 10*G(:,j).*A;
        dB = p*A - m*B.*(1-A);
    end

end