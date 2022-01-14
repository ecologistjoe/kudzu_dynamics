function O = bestPulsesInvasion(N)

[b1 b2] = meshgrid([0 .25 1 2 5], [0 .25 1 2 5]);
B0 = [b1(:) b2(:)];
B0 = B0(mod(N-1,25)+1,:);
k = 0;

switch ceil(N/25)
case 1                  %A self contained field. Only interested in having half clear.
    S = [-1  1
          1 -1]/2;
case 2                  %Two plots with a small connection.
    S = [-1  1
          1 -1]/6;
case 3                  %A plot connected to a large, agressive field by a 'bridge' or 'corridor' (plot 2)
    S = [-1  1
          2  0]/6;
end 


V{1} = nchoosek(1:240, 1);
V{2} = nchoosek(1:240, 2);

for Method = [0 3]
for m = [1 2 3 5]
for p = [1 2 3 5]
    if(p<m), continue; end;
for days=1:20
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
            [x, J] = chooseDaysGA(days, 2, 30, 4, m, p, B0, S, Method,0);
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
    R.time = toc(start);
    O(k) = R;
end
   %save(sprintf('output_3_%d.dat',N), 'O');
   save(sprintf('/data/scratch/mhughe13/output_invasion_%d.dat',N), 'O');
end
end
end