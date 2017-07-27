%% Load data

% clear all
% close all
% 
% [nodePot2,edgePot2,edgeStruct,nRows,nCols,maxVal] = middlebury_loadDataUGM('tsukuba');
% [nNodes,nStates] = size(nodePot2);
% 
% w = nRows;
% h = nCols;
% [edges_out,edges_in] = get_directed_struct(edgeStruct);
% 
% 
% 
% nEdges = size(edgePot2,3);
% d = nStates;
% 
% 
% nodePot = -log(nodePot2);
% edgePot = -log(edgePot2);
% 
% one = ones(d,1);
% C = repmat(diag(one),1,d);
% D = kron(eye(d),one');
% 
% A = C'*C + D'*D;

% S = d^2*nEdges;
% Q = sparse(S, S); %// preallocates with zeros
% for i = 0:d^2:S-1
%     Q(i+(1:d^2),i+(1:d^2)) = A;
% end

MSparse = sparse(A);
% replicate: please don't test with n=220512 directly
Mblk = repmat({MSparse}, 1, nEdges);
% convert do block diagonal
Q = blkdiag(Mblk{:});

lb = zeros(d*d*nEdges,1);
ub = ones(d*d*nEdges,1);
x0 = lb;




theta_nodes = nodePot';
theta_edges = reshape(edgePot,d*d,nEdges);

x_nodes = zeros(d,nNodes);
x_edges = zeros(d*d,nEdges);

y_edges = zeros(d,nEdges);
z_edges = zeros(d,nEdges);

%y_nodes = zeros(d,nNodes);
%z_nodes = zeros(d,nNodes);

avg_out = zeros(d,nNodes);
avg_in = zeros(d,nNodes);

delta_out = zeros(d,nEdges);
delta_in = zeros(d,nEdges);

lambda_nodes = zeros(d,nNodes);
lambda_edges = zeros(d,nEdges);

mu_nodes = zeros(d,nNodes);
mu_edges = zeros(d,nEdges);

B = zeros(d*d,nEdges);



%% Decode by DD

dual_current    = 0;
dual_bestsofar  = 0;
primal_current  = 0;

dual_best = 484.7098;

new_lambda_nodes = lambda_nodes;
new_mu_nodes = mu_nodes;
new_lambda_edges = lambda_edges;
new_mu_edges = mu_edges;

rho = 1;

dual = [];
primal = [];


% i = 0;
for i= 1:1000
    fprintf('\n%d.\n',i);
    %tStart = tic; 
    
    lambda_nodes = new_lambda_nodes;
    mu_nodes = new_mu_nodes;
    lambda_edges = new_lambda_edges;
    mu_edges = new_mu_edges;
    
    
    %%
    %% 1. Solve the node slaves
    %%
    
    
    fprintf('Solve node slaves... \n')
    for p = 1:nNodes
        %%% Compute y_p and z_p
        edges1 = edges_out{p};
        edges2 = edges_in{p};
        y_p = -sum(y_edges(:,edges1),2);        
        z_p = -sum(z_edges(:,edges2),2);
        
        a = 2/(rho*(length(edges1)+length(edges2)))*(theta_nodes(:,p) + lambda_nodes(:,p) + mu_nodes(:,p) + rho*y_p + rho*z_p);
        k = (sum(a)+2)/d;
        
        x_nodes(:,p) = 0.5*(k-a);
    end    
    
    
        
    %% Solve the edge slaves
    %tic   
       
%     fprintf('Optimize... \n')
%     for e = 1:nEdges
%         if(rem(e,100) == 0)
%             fprintf('%d/%d\n',e,nEdges);
%         end
%         p = edgeStruct.edgeEnds(e,1);
%         q = edgeStruct.edgeEnds(e,2);
%         
%         b = (1/rho)*(theta_edges(:,e)) + (1/rho)*(D'*lambda_edges(:,e)) + (1/rho)*(C'*mu_edges(:,e))...
%              - D'*x_nodes(:,p) - C'*x_nodes(:,q);
%         
% %         options = optimoptions(@quadprog,'Algorithm','active-set');
% %         x = quadprog(A,b,[],[],[],[],lb,ub,ub,options); 
%         
%         cvx_begin quiet,
%             variable x(d*d,1)
%             x >= 0
%             minimize(0.5*x'*A*x + b'*x)
%         cvx_end
%         x_edges(:,e) =  x;        
%     end
    
    
    %% Compute B = [b_1 ... b_pq ... b_nEdges]
    fprintf('Compute B... \n')
    for e = 1:nEdges
        p = edgeStruct.edgeEnds(e,1);
        q = edgeStruct.edgeEnds(e,2);
        
        B(:,e) = (1/rho)*(theta_edges(:,e)) + (1/rho)*(D'*lambda_edges(:,e)) + (1/rho)*(C'*mu_edges(:,e))...
            - D'*x_nodes(:,p) - C'*x_nodes(:,q);
              
    end
    
    
    c = reshape(B,d*d*nEdges,1);
    
    
    fprintf('Optimize... \n')
    
    options = optimoptions(@quadprog,'Algorithm','interior-point-convex');
    x = quadprog(Q,c,[],[],[],[],lb,ub,[],options);
    
%     cvx_begin quiet
%         variable X(d*d,nEdges)
%         minimize(0.5*sum(sum_square(D*X)) + 0.5*sum(sum_square(C*X)) + sum(sum(X.*B)));
%         subject to
%             X >= 0;
%     cvx_end       

    x_edges =  reshape(x,d*d,nEdges);  
    
    clear x;
    
    
    fprintf('Compute dual gap... \n')
    %% Compute necessary quantities
    
    %%% Compute y_pq = D*x_pq and z_pq = C*x_pq
    M =  reshape(x_edges,[d,d,nEdges]);    
    
    y_edges = squeeze(sum(M,1));
    z_edges = squeeze(sum(M,2));   
    
    x_nodes_vec = reshape(x_nodes,d*nNodes,1);
    y_edges_vec = reshape(y_edges,d*nEdges,1);
    z_edges_vec = reshape(z_edges,d*nEdges,1);
    
              
    dual_current = x_nodes_vec'*reshape(theta_nodes + lambda_nodes + mu_nodes,d*nNodes,1) + reshape(x_edges,d*d*nEdges,1)'*reshape(theta_edges + D_lambda_edges + C_mu_edges,d*d*nEdges,1);
    
    [~,ind] = max(x_nodes,[],1);
    
    primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    
    dual(i) = dual_current;
    primal(i) = primal_current;
    
    fprintf('Dual: %d, Primal: %d',dual_current,primal_current);
    
    if(dual_best - dual_current < 1)
        break;
    end
    
    
    %% %% 
    %% %% Update the parameters    
    %% %% 
    
    %%% Compute the averages
    %get_avg(avg_out,avg_in,x_nodes,y_edges,z_edges,edges_out,edges_in);
    
    for p = 1:size(x_nodes,2)
        edges1 = edges_out{p};
        s1 = sum(y_edges(:,edges1),2);
        avg_out(:,p) = ( 1/(length(edges1)+1) )* ( x_nodes(:,p) + s1 );

        edges2 = edges_in{p};
        s2 = sum(z_edges(:,edges2),2);
        avg_in(:,p) = ( 1/(length(edges2)+1) )* ( x_nodes(:,p) + s2 );
    end
    
    
    new_lambda_nodes    = lambda_nodes + rho*(x_nodes - avg_out);
    new_mu_nodes        = mu_nodes + rho*(x_nodes - avg_in);
    
    %toc
    
    for e = 1:nEdges
        p = edgeStruct.edgeEnds(e,1);
        q = edgeStruct.edgeEnds(e,2);
        delta_out(:,e) = y_edges(:,e) - avg_out(:,p);
        delta_in(:,e) = z_edges(:,e) - avg_in(:,q);
    end
    
    %%% Update theta_pq = theta_pq + alpha*D'*delta_out + alpha*C'*delta_in
    %theta_edges = theta_edges + alpha*delta_out(kron((1:d)',ones(d,1)),:) + alpha*repmat(delta_in,d,1);
    new_lambda_edges    = lambda_edges + rho*delta_out;
    new_mu_edges        = mu_edges + rho*delta_in;
    
    
end


Energy_DD = -UGM_LogConfigurationPotential(ind,nodePot2,edgePot2,edgeStruct.edgeEnds);

figure;imagesc(reshape(ind,nRows,nCols)');title(sprintf('DD (Energy = %f)',Energy_DD));colormap gray;pause(1);
print('-depsc','results');


% save('dual','dual');
% save('primal','primal');

% i_primal = 1:length(primal);
% i_dual = 1:length(dual);
% 
% 
% 
% figure;plot(i_primal, primal,i_dual, dual);
% xlim([0,4201]);
% xlabel('Iterations')
% ylabel('Energy');
% hleg1 = legend('Primal','Dual');
% title('Dual decomposition');
% print('-depsc','energy');
