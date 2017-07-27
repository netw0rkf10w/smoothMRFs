function [edges_out,edges_in] = get_directed_struct(edgeStruct)

    nNodes = edgeStruct.nNodes;
    
    edges_out    = cell(1,nNodes);
    edges_in     = cell(1,nNodes);    
    
    for p = 1:nNodes
        %%% Check all the edges associated to the node p
        edges = UGM_getEdges(p,edgeStruct);
        edges_out{p} = [];
        edges_in{p} = [];
        %%% For each edge, check if it is "out" or "in"
        for i=1:length(edges)
            %%% Get the neighbor q
            nodes = edgeStruct.edgeEnds(edges(i),:);
            q = nodes(nodes ~= p);
            
            if(q > p)
                edges_out{p} = [edges_out{p} edges(i)];
            else
                edges_in{p} = [edges_in{p} edges(i)];
            end
        end      
    end
end