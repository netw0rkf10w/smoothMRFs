function get_avg(avg_out,avg_in,x_nodes,y_edges,z_edges,edges_out,edges_in)

for p = 1:size(x_nodes,2)
    edges1 = edges_out{p};
    s1 = sum(y_edges(:,edges1),2);
    avg_out(:,p) = ( 1/(length(edges1)+1) )* ( x_nodes(:,p) + s1 );
    
    edges2 = edges_in{p};
    s2 = sum(z_edges(:,edges2),2);
    avg_in(:,p) = ( 1/(length(edges2)+1) )* ( x_nodes(:,p) + s2 );
end

end
