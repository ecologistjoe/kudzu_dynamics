m = 2;
p = 2;
B0 = .5;

S = struct();
for k=1:4
    if(k <= 5)
        V = nchoosek(1:30, k);
    else
        V = load(sprintf('nc%d.mat', k));
    end
    
    J = scorePulse(V, 30, 1, m, p, B0, 0, 0,1);
    Jmin = min(J);
    Best = V(find(J==Jmin),:);
    
    S(k).J = Jmin;
    S(k).V = Best;
end
   