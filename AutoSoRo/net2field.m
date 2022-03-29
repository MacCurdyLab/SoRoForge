function [Q] = net2field(Net,Params,varargin)

%Net is a struct with the following fields
%n = a Nx1 vector of node names, where N is the number of nodes
%b = a Nx1 vector of node biases, where N is the number of nodes
%a = a Nx1 vector of node aggregation type, where N is the number of nodes

%s = a Mx1 vector of startnodes, where M is the number of links
%t = a Mx1 vector of termination nodes, where M is the number of links
%w = a Mx1 vector of link weights, where M is the number of links

%build the directed graph from the matrix representation
G = digraph(Net.s,Net.t,Net.w);

if isempty(varargin)
%build the coorndinate field
[X,Y,Z] = meshgrid(Params.Extents(1,1):Params.dx:Params.Extents(1,2),...
    Params.Extents(2,1):Params.dx:Params.Extents(2,2),...
    Params.Extents(3,1):Params.dx:Params.Extents(3,2));
else
   C = varargin{1}./10;
   X = C(:,1); Y = C(:,2); Z = C(:,3);
end

%produce a set of cylindrical coordinates
[T,R,Z] = cart2pol(X,Y,Z);

%Initialize the first five (default) nodes and their biases
nodevals = cell(length(Net.n),1);
nodevals{1} = X + Net.b(1); 
nodevals{2} = Y + Net.b(2); 
nodevals{3} = Z + Net.b(3);
nodevals{4} = T + Net.b(4); 
nodevals{5} = R + Net.b(5);

%Initialize the rest of the node values
for i = 6:length(Net.n)
    nodevals{i} = Net.b(i);
end

%extract the cell array of node names, to be whittled away
nodes = Net.n;

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
            operation = string(nodes(non_startpts(i)));
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

function Q = computeNode(name,inputs,weights,offset,aggtype)
    
    %initialize the aggregated node with the first input
    aggregated = inputs{1}*weights(1);
    
    %check if we need to aggregate
    if length(inputs)>1
       if contains(name,'*') %aggragate by multiplication
           name = strrep(name,'*','');
           for i = 2:length(inputs)
                aggregated = aggregated .*inputs{i}*weights(i);
           end
       else                      %aggragate by addition
           for i = 2:length(inputs)
                aggregated = aggregated + inputs{i}*weights(i);
           end
       end
    end
        
    if strcmp(name,'i') || contains(name,'Identity')
        Q = aggregated               + offset;
    elseif strcmp(name,'a') || contains(name,'Atan')
        Q = atan(aggregated)          + offset;
    elseif strcmp(name,'s')  || contains(name,'Sine')               
        Q = sin(aggregated)          + offset;
    elseif strcmp(name,'c')  || contains(name,'Cosine')               
        Q = cos(aggregated)          + offset;
    elseif strcmp(name,'f') || contains(name,'Fourth')
        Q = aggregated.^4            + offset;
    elseif strcmp(name,'q') || contains(name,'Squared')
        Q = aggregated.^2            + offset;
    elseif strcmp(name,'n') || contains(name,'Asin')
        Q = asin(aggregated)            + offset;
    elseif strcmp(name,'o') || contains(name,'Acos')
        Q = acos(aggregated)            + offset;
    elseif strcmp(name,'p') || contains(name,'Inv')
        Q = (aggregated).^(-1)            + offset;
   elseif strcmp(name,'r') || contains(name,'Sqrt')
        Q = (aggregated).^(0.5)            + offset;
    end
end