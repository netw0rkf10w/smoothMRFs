function [func, grad_lambda, grad_nu, u_nodes, u_edges] =  compute(lambda, nu, mu, theta_nodes, theta_edges, edgeStruct, edges_out, edges_in)
    
    %% Compute u(x)
    
    %% Compute the gradient by solving the node and edge subproblems
    %fprintf('Solve node slaves... \n')
    
    [d,nNodes] = size(theta_nodes);
    nEdges = size(theta_edges,2);   
   
    lambda_nodes = zeros(d,nNodes);
    nu_nodes = zeros(d,nNodes);
    parfor p = 1:nNodes
        %%% Compute v_p and w_p
        edges1 = edges_out{p};
        edges2 = edges_in{p};
        lambda_nodes(:,p) = -sum(lambda(:,edges1),2);        
        nu_nodes(:,p) = -sum(nu(:,edges2),2);        
    end   
    
    a = -1/mu*(theta_nodes + lambda_nodes + nu_nodes);   
    expAMinusMax = exp(bsxfun(@minus,a,max(a)));
    sumExpAMinusMax = sum(expAMinusMax);
    u_nodes = bsxfun(@rdivide,expAMinusMax,sumExpAMinusMax);  
    
    
    %% Solve the edge slaves
    %fprintf('Solve edge slaves... \n')
    %%% Compute D'*lambda_pq + C'*nu_pq
    D_lambda = lambda(kron((1:d)',ones(d,1)),:);
    C_nu = repmat(nu,d,1); 
    
    b = -1/mu*(theta_edges + D_lambda + C_nu);
    expBMinusMax = exp(bsxfun(@minus,b,max(b)));
    sumExpBMinusMax = sum(expBMinusMax);
    
    %%% Divide each column by its sum
    u_edges = bsxfun(@rdivide,expBMinusMax,sumExpBMinusMax);    
    
    
    %% Compute the value of f_mu(x) = f_mu(lambda,nu)
    func = sum(max(a)) + sum(log(sumExpAMinusMax)) + sum(max(b)) + sum(log(sumExpBMinusMax)) - (nNodes + 2*nEdges)*log(d);
    func = mu*func;
    
    %% Compute the gradient g(lambda,nu)
    
    %% First compute necessary quantities
    
    %%% Compute v_pq = D*u_pq and w_pq = C*u_pq
    M =  reshape(u_edges,[d,d,nEdges]);
    v_edges = squeeze(sum(M,1));
    w_edges = squeeze(sum(M,2));   
    
    
    parfor e = 1:nEdges
        p = edgeStruct.edgeEnds(e,1);
        q = edgeStruct.edgeEnds(e,2);
        grad_lambda(:,e) = -v_edges(:,e) + u_nodes(:,p);
        grad_nu(:,e) = -w_edges(:,e) + u_nodes(:,q);
    end    
    
    
    %% Compute the value of f_mu(x) = f_mu(lambda,nu)
%     func =  reshape(lambda,d*nEdges,1)'*reshape(grad_lambda,d*nEdges,1)...
%         + reshape(nu,d*nEdges,1)'*reshape(grad_nu,d*nEdges,1)...
%         - reshape(theta_nodes,d*nNodes,1)'*reshape(u_nodes,d*nNodes,1)...
%         - reshape(theta_edges,d*d*nEdges,1)'*reshape(u_edges,d*d*nEdges,1)...
%         - mu*sum(sum(u_nodes.*log(u_nodes+(u_nodes==0)))) + mu*sum(sum(u_edges.*log(u_edges+(u_edges==0))))...
%         - mu*(nNodes+2*nEdges)*log(d); 
    
    
end

