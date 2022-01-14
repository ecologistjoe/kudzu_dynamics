function r= numfetch(conn, q)
    r = fetch(conn, q);
    if(isempty(r)), return; end;
    fn = fieldnames(r);
    for i=1:length(fn);
        f = cell2mat(fn(i));
        if(iscell(r.(f)))
            try
                r.(f) = cellfun(@str2num, r.(f));
            catch
                r.(f) = cellfun(@str2num, r.(f), 'UniformOutput', false);
            end
        end
    end
end


