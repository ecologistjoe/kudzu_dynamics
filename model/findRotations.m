function Z = findRotations(L)

R =[2   3   4   5   6   7
    3   4   5   6   7   2
    4   5   6   7   2   3
    5   6   7   2   3   4
    6   7   2   3   4   5
    7   2   3   4   5   6
    7   6   5   4   3   2
    6   5   4   3   2   7
    5   4   3   2   7   6
    4   3   2   7   6   5
    3   2   7   6   5   4
    2   7   6   5   4   3]-1;

N = length(L);
p = N.^(0:5);
pr = repmat(p, 12, 1);

Si = [];
S = [];
z = 0;

for a = 1:N
    L(a)
    for b = 1:N
        for c = 1:N
            if(sum(L([a b c])) > 6) continue; end;
            for d = 1:N
                for e = 1:N
                    for f = 1:N                       
                        %Form the permutation and all of its rotations
                        T = [a b c d e f];
                        if(sum(L(T)) == 6)
                            %the max works as the index because we add higher
                            %values first
                            idx = max(sum(T(R).*pr, 2));

                            if(all(idx ~= Si))
                                z = z+1;
                                S(z,:) = L(T);      %add the values of this permutation
                                Si(z) = idx;
                            end  
                        end
                    end
                end
            end
        end
    end
end
Z = [];
for a = L
    a = repmat(a, length(S), 1);
    Z = [Z; a S];
end     