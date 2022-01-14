function [J varargout] = chooseGoatDays(x,sendToWorkspace,w,IC)
%% Determine the Form of inputs and Create the Control Regime
    if(nargin<2)
        sendToWorkspace = 0;
    end

    [rows cols] = size(x);
    
    if(cols ==300)
        x = x(:,floor(1:.5:300.5));
    end
    
    G1 = x(:,1:200)*10;
    G2 = x(:,201:400)*10;
    G3 = x(:,401:600)*10;

    
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
    
    A1 = ones(rows,1) * IC(1);
    B1 = ones(rows,1) * IC(2);
    
    A2 = ones(rows,1) * IC(3);
    B2 = ones(rows,1) * IC(4);

    A3 = ones(rows,1) * IC(5);
    B3 = ones(rows,1) * IC(6);
    
    m = 3;
    p = 3;
    N = 200;
    n = 20;
    h = 1/(n*N);
    t = h;
    h2 = h/2;
    h6 = h/6;
    A1s = zeros(rows, N);
    B1s = zeros(rows, N);
    A2s = zeros(rows, N);
    B2s = zeros(rows, N);
    A3s = zeros(rows, N);
    B3s = zeros(rows, N);
    
%% Solve ODE System
    for j = 1:N
        A1s(:,j)=A1;
        B1s(:,j)=B1;
        A2s(:,j)=A2;
        B2s(:,j)=B2;
        A3s(:,j)=A3;
        B3s(:,j)=B3;
        
        for i = 1:n
            [k1A1 k1B1 k1A2 k1B2 k1A3 k1B3] = growAB(A1,          B1,          A2,          B2,          A3,          B3,  t);
            [k2A1 k2B1 k2A2 k2B2 k2A3 k2B3] = growAB(A1+h2.*k1A1, B1+h2.*k1B1, A2+h2.*k1A2, B2+h2.*k1B2, A3+h2.*k1A3, B3+h2.*k1B3,  t);
            [k3A1 k3B1 k3A2 k3B2 k3A3 k3B3] = growAB(A1+h2.*k2A1, B1+h2.*k2B1, A2+h2.*k2A2, B2+h2.*k2B2, A3+h2.*k2A3, B3+h2.*k2B3,  t);
            [k4A1 k4B1 k4A2 k4B2 k4A3 k4B3] = growAB(A1+h .*k3A1, B1+h .*k3B1, A2+h .*k3A2, B2+h .*k3B2, A3+h .*k3A3, B3+h .*k3B3,  t);
            A1 = A1 + h6 * (k1A1 + 2*k2A1 + 2*k3A1 + k4A1);
            B1 = B1 + h6 * (k1B1 + 2*k2B1 + 2*k3B1 + k4B1);
            A2 = A2 + h6 * (k1A2 + 2*k2A2 + 2*k3A2 + k4A2);
            B2 = B2 + h6 * (k1B2 + 2*k2B2 + 2*k3B2 + k4B2);
            A3 = A3 + h6 * (k1A3 + 2*k2A3 + 2*k3A3 + k4A3);
            B3 = B3 + h6 * (k1B3 + 2*k2B3 + 2*k3B3 + k4B3);
            t = t + h;
        end
    end
    
%% Find Value of the Objective Function     
    if(nargin<3)
        w = [3000 200 100 1];
    end
    
    o = zeros(rows,1);

    %Projected Amount of Future Aboveground Growth
    R1 = 2 * log(1+1./(5*B1.^2)) ./ (m*(sqrt(5)-1)); % From Maple, see notes.
    R2 = 2 * log(1+1./(5*B2.^2)) ./ (m*(sqrt(5)-1)); % From Maple, see notes.
    R3 = 2 * log(1+1./(5*B3.^2)) ./ (m*(sqrt(5)-1)); % From Maple, see notes.
    
    J = w(1)*(3-min(R1,2) - sum(A1s<.1,2)/N) + ...
        w(2)*sum(diff([o G1 o],1,2)~=0,2) + ...
        w(3)*sum(G1>0,2) + ...
        w(4)*sum(max(G1-100*A1s,0),2) +...
        ...
        w(2)*sum(diff([o G2 o],1,2)~=0,2) + ...
        w(3)*sum(G2>0,2) + ...
        w(4)*sum(max(G2-100*A2s,0),2) +...
        ...
        w(2)*sum(diff([o G3 o],1,2)~=0,2) + ...
        w(3)*sum(G3>0,2) + ...
        w(4)*sum(max(G3-100*A3s,0),2);
    
    J(J<0) = 1e10;
    J(isnan(J)) = 1e10;

%% Place some values in the Workspace if Requested
    if(sendToWorkspace)
        varargout(1) = {[A1s' A2s' A3s']};
        varargout(2) = {[B1s' B2s' B3s']};
        varargout(3) = {[G1' G2' G3']};
    end
    
%% Define the Growth Function
    function [dA1 dB1 dA2 dB2 dA3 dB3] = growAB(A1, B1, A2, B2, A3, B3, t)
        dA1 = m*B1.*(1-A1) - 10*G1(:,j).*A1 + 0.5*( -A1 +   A2);
        dA2 = m*B2.*(1-A2) - 10*G2(:,j).*A2 + 0.5*(  A1 - 2*A2 + A3);
        dA3 = m*B3.*(1-A3) - 10*G3(:,j).*A3 + 0.5*(         A2 - A3);
        
        dB1 = m*A1 - m*B1.*(1-A1);
        dB2 = m*A2 - m*B2.*(1-A2);
        dB3 = m*A3 - m*B3.*(1-A3);
    end

end