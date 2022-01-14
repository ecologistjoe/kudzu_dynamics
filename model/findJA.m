conn = database('kudzu','root','','com.mysql.jdbc.Driver',...
                'jdbc:mysql://localhost/kudzu');

S=[-6   1   1   1   1   1   1
    1  -3   1   0   0   0   1
    1   1  -3   1   0   0   0
    1   0   1  -3   1   0   0
    1   0   0   1  -3   1   0
    1   0   0   0   1  -3   1
    1   1   0   0   0   1  -3]/12;

load rotations;
B0 = rotations;
clear rotations;

for pulses=[1 2 3 4 6 8 10 15 20]
    pulses
    for method=[0 3]
        q = sprintf(['SELECT id, m, p, b0, x FROM bestpulses '...
                    ' WHERE run="INV_2D" AND pulses=%d AND method=%d'], pulses, method);
        R = numfetch(conn, q);

        if(pulses>1); R.x = cell2mat(R.x); end

        if(method == 0)
            meth = 3;
        else
            meth = 0;
        end
           meth=method;
        J = scorePulse(R.x, 30, 4, R.m, R.p, cell2mat(R.b0), S, meth,0);


        v = sprintf('(%d, %.32f),\n', [R.id J]');
        v = v(1:end-2);


        r = exec(conn, [sprintf('INSERT IGNORE bestpulses (id, score%d) VALUES ', meth) v ...
                        sprintf(' ON DUPLICATE KEY UPDATE score%d = VALUES(score%d)', meth,meth)]);
    end
end

