function [Nodes] = makeNodeNamesUnique(Nodes)

for i = 1:length(Nodes)
   if contains('0123456789',Nodes{i}(end))
       Nodes{i}(end) = [];
   end
end


cache = {};
count = [];
for i = 1:length(Nodes)
    [Li, Lo] = ismember(Nodes{i},cache);
    if ~Li
        cache{end+1}=Nodes{i};
        count(end+1)=1;
    else
        count(Lo) = count(Lo)+1;
        newname = [Nodes{i} num2str(count(Lo))];
        Nodes{i} = newname;
    end
end

end