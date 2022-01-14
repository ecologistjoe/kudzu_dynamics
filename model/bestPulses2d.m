function O = bestPulses2d(N)

load rotations;
B0 = rotations(N,:);

%Define the Spread Matrix.  Cell 1 is central, cells 2-7 surround it.
S=[-6   1   1   1   1   1   1
    1  -3   1   0   0   0   1
    1   1  -3   1   0   0   0
    1   0   1  -3   1   0   0
    1   0   0   1  -3   1   0
    1   0   0   0   1  -3   1
    1   1   0   0   0   1  -3]/12;

k = 0;


V{1} = nchoosek(1:840, 1);
V{2} = nchoosek(1:840, 2);

%Do 20 combinations of Methods, m, and p.
for Method = [0 3]
for m = [1 2 3 5]
for p = [1 2 3 5]
    if(p<m), continue; end;
for days=[1 2 3 4 6 8 10 15 20]
    start = tic;
    k = k+1;
    R = struct();
    if(days<=2)
        v = V{days};
        J = scorePulse(v, 30, 4, m, p, B0, S, Method, 0);
        R.J = min(J);
        R.x = v(J==R.J,:);
        R.x = R.x(1,:);
    else
        R.J = inf;
        for i=1:3
            [x, J] = chooseDaysGA(days, 7, 30, 4, m, p, B0, S, Method,0);
            if(J < R.J)
                R.J=J(1);
                R.x=x(1,:);
            end
        end
    end
    
    R.J0 = scorePulse(R.x, 30, 4, m, p, B0, S, Method, 0);
    R.B0 = B0;
    R.days = days;
    R.m = m;
    R.p = p;
    R.Method = Method;
    R.time = toc(start);
    O(k) = R;
    
end
   %save(sprintf('output_3_%d.dat',N), 'O');
   save(sprintf('/data/scratch/mhughe13/output_3_%d.dat',N), 'O');
end
end
end
