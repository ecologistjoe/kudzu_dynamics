function O = bestPulsesBaseLineRandom()


k = 0;

V{1} = nchoosek(1:120, 1);
V{2} = nchoosek(1:120, 2);

for B0 = [.25 1 2 5]
for Method = [0 3]
for m = [1 2 3 5]
for p = [1 2 3 5]
    if(p<m), continue; end;
disp([B0 Method m p]);
    
for days=1:20
    start = tic;
    k = k+1;
    R = struct();
    if(days<3)
        v = V{days};
        J = scorePulse(v, 30, 4, m, p, B0, 0, Method, 0);
    else
        v = [];
        n = 5000;
        while (size(v,1) < n)
            r= randi([1 120]-1, 2*n, days)*2*n + repmat((1:(2*n))', 1, days);
            
            z = zeros(2*n, 120);
            z(r(:)) = 1;
            z(sum(z, 2) < days,:) = [];
            v = [v;z];
        end
        v = logical(v(1:n,:));
        J = scorePulse(v, 30, 4, m, p, B0, 0, Method, 0);
    end
    
    R.B0 = B0;
    R.days = days;
    R.m = m;
    R.p = p;
    R.Method = Method;
    R.time = toc(start);
    
    R.minJ = min(J);
    R.maxJ = max(J);
    R.meanJ = mean(J);
    R.medJ = median(J);
    R.stdJ = std(J);
    
    O(k) = R;
end
   save('output_rnd.dat', 'O');
   %save(sprintf('/data/scratch/mhughe13/output_seasonal_%d.dat',N), 'O');
end
end
end
end