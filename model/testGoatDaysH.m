function testGoatDaysH()

G = zeros(1,250);
G([60 61 62 80 89 239]) = 1;
T = length(G);
Bs = zeros(10,T);
As = zeros(10,T);
    
for i = 1:10
    
%% Set up Problem
    A = 0.05;
    B = .5;
    z = 1;
    h = i/10;
    t = h;
    h2 = h/2;
    h6 = h/6;
    j=0;
 
%% Solve ODE System
    while (t<=T)
        j=j+1;
        Bs(i, ceil(t))=B;
        As(i, ceil(t))=A;

        [k1A k1B] = growAB(A, B, t);
        [k2A k2B] = growAB(A+h2*k1A, B+h2*k1B, t);
        [k3A k3B] = growAB(A+h2*k2A, B+h2*k2B, t);
        [k4A k4B] = growAB(A+h*k3A, B+h*k3B, t);
        A = A + h6 * (k1A + 2*k2A + 2*k3A + k4A);
        B = B + h6 * (k1B + 2*k2B + 2*k3B + k4B);
        t = t + h;
    end
end

%% Place some values in the Workspace if Requested
        assignin('base', 'B', Bs);
        assignin('base', 'A', As);
        assignin('base', 'G', G);
    
    
%% Define the Growth Function
    function [dA dB] = growAB(A, B, t)
        dA = 0.05*B*(1-A) - G(ceil(t))*A;
        dB = 0.05*A*(1-B) - 0.05*B*(1-A);
    end
end
