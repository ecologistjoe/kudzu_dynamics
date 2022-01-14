function O = bestPulsesSeasonal(N)

B0 = [.25 1 2 5];
B0 = B0(N);
k = 0;

V{1} = nchoosek(1:120, 1);
V{2} = nchoosek(1:120, 2);
V{3} = nchoosek(1:120, 3);

for Method = [0]
for m = [2]
for p = [2]
    if(p<m), continue; end;
for days=1:4
    start = tic;
    k = k+1;
    [k p m days]
    R = struct();
    if(days<=3)
        v = V{days};
        J = scorePulse(v, 30, 4, m, p, B0, 0, Method, 0);
        R.J = min(J);
        R.x = v(J==R.J,:);
        R.x = R.x(1,:);
    else
        R.J = inf;
        for i=1:3
            [x, J] = chooseDaysGA(days, 1, 30, 4, m, p, B0, 0, Method,0);
            if(J < R.J)
                R.J=J(1);
                R.x=x(1,:);
            end
        end
    end
    
    R.B0 = B0;
    R.days = days;
    R.m = m;
    R.p = p;
    R.Method = Method;
    R.time = toc(start)
    O(k) = R;
end
   %save(sprintf('output_3_%d.dat',N), 'O');
   save(sprintf('/data/scratch/mhughe13/output_seasonal_%d.dat',N), 'O');
end
end
end