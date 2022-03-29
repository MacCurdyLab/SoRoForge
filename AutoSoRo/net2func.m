function strout = net2func(Net)

%Net is a struct with the following fields
%n = a Nx1 vector of node names, where N is the number of nodes
%b = a Nx1 vector of node biases, where N is the number of nodes
%a = a Nx1 vector of node aggregation type, where N is the number of nodes

%s = a Mx1 vector of startnodes, where M is the number of links
%t = a Mx1 vector of termination nodes, where M is the number of links
%w = a Mx1 vector of link weights, where M is the number of links

%build the directed graph from the matrix representation
G = digraph(Net.s,Net.t,Net.w);

%extract the cell array of node names, to be whittled away
nodes = Net.n;

%Initialize the first five (default) nodes and their biases
nodevals = cell(length(Net.n),1);
nodevals{1} = ['(x+' num2str(Net.b(1)) ')']; 
nodevals{2} = ['(y+' num2str(Net.b(2)) ')']; 
nodevals{3} = ['(z+' num2str(Net.b(3)) ')'];
nodevals{4} = ['(t+' num2str(Net.b(4)) ')']; 
nodevals{5} = ['(r+' num2str(Net.b(5)) ')'];

%Initialize the rest of the node values
for i = 6:length(Net.n)
    nodevals{i} = num2str(Net.b(i));
end

%intialize string
strout = [];

while numedges(G)>0
    
    %remove islands from graph
    islands = findIslands(G);
    G = rmnode(G,islands);
    nodes(islands) = [];
    nodevals(islands) = [];

    %find startpts
    startpts = findStarts(G);
    
    %all the other nodes are not startnodes
    non_startpts = 1:size(G.Nodes,1); 
    non_startpts(startpts) = [];
    
    i = 0;
    while i<length(non_startpts)
        i=i+1;
        [ein, nid] = inedges(G,non_startpts(i));
        weights = G.Edges.Weight(ein);
        offset = nodevals{non_startpts(i)};       
        if all(ismember(nid,startpts))

            %compute the new value of this node
            operation = nodes{non_startpts(i)};
            inputs = nodevals(nid);
            Q = computeNode(operation,inputs,weights,offset);

            %Assign this new value to the node
            nodevals{non_startpts(i)} = Q;

            %remove the edges that were inputs
            G = rmedge(G,nid,ones(length(nid),1)*non_startpts(i));
            i = length(non_startpts);              
        end      
    end
    
end
   strout = string(str2sym(Q)); 
end

function islands = findIslands(G)
    islands = [];
    for i = 1:size(G.Nodes,1)
       if isempty(outedges(G,i)) && isempty(inedges(G,i))
           islands = [islands i];
       end
    end
end

function startpts = findStarts(G)
    startpts = [];
    for i = 1:size(G.Nodes,1)
       if ~isempty(outedges(G,i)) && isempty(inedges(G,i))
           startpts = [startpts i];
       end
    end
end

function Q = computeNode(name,inputs,weights,offset)

    %initialize the aggregated node with the first input
    Q = [inputs{1} '*' num2str(weights(1))] ;
    
    %check if we need to aggregate
    if length(inputs)>1
       if contains(name,'*') %aggragate by multiplication
           name = strrep(name,'*','');
           for i = 2:length(inputs)
                Q = [Q '*' inputs{i} '*' num2str(weights(i))] ;
           end
       else                      %aggragate by addition
           for i = 2:length(inputs)
                Q = [Q '+' inputs{i} '*' num2str(weights(i))] ;
           end
       end
    end
        
    if strcmp(name,'i') || contains(name,'Identity')
        Q = [Q '+' offset];
    elseif strcmp(name,'a') || contains(name,'Atan')
        Q = ['(atan(' Q ')+' offset ')'];
    elseif strcmp(name,'s')  || contains(name,'Sine')               
        Q = ['(sin(' Q ')+' offset ')'];
    elseif strcmp(name,'c')  || contains(name,'Cosine')               
        Q = ['(cos(' Q ')+' offset ')'];
    elseif strcmp(name,'f') || contains(name,'Fourth')
        Q = ['((' Q ')^4 +' offset ')'];
    elseif strcmp(name,'q') || contains(name,'Squared')
        Q = ['((' Q ')^2 +' offset ')'];
    elseif strcmp(name,'n') || contains(name,'Asin')
        Q = ['(asin(' Q ')+' offset ')'];
    elseif strcmp(name,'o') || contains(name,'Acos')
        Q = ['(acos(' Q ')+' offset ')'];
    elseif strcmp(name,'p') || contains(name,'Inv')
        Q = ['((' Q '^(-1))+' offset ')'];
    elseif strcmp(name,'r') || contains(name,'Sqrt')
        Q = ['(sqrt(' Q ')+' offset ')'];
    end
end


