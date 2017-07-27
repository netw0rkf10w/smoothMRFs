
%% Load data

clear all
close all

[nodePot2,edgePot2,edgeStruct,nRows,nCols,maxVal] = middlebury_loadDataUGM('tsukuba');
[nNodes,nStates] = size(nodePot2);

w = nRows;
h = nCols;
[edges_out,edges_in] = get_directed_struct(edgeStruct);



nEdges = size(edgePot2,3);
d = nStates;


nodePot = -log(nodePot2);
edgePot = -log(edgePot2);




theta_nodes = nodePot';
theta_edges = reshape(edgePot,d*d,nEdges);

x_nodes = zeros(d,nNodes);
x_edges = zeros(d*d,nEdges);

y_edges = zeros(d,nEdges);
z_edges = zeros(d,nEdges);

avg_out = zeros(d,nNodes);
avg_in = zeros(d,nNodes);

delta_out = zeros(d,nEdges);
delta_in = zeros(d,nEdges);

lambda_nodes = zeros(d,nNodes);
lambda_edges = zeros(d,nEdges);

mu_nodes = zeros(d,nNodes);
mu_edges = zeros(d,nEdges);


%% Decode by DD

dual_current    = 0;
dual_bestsofar  = 0;
primal_current  = 0;

dual_best = 484.7098;

new_lambda_nodes = lambda_nodes;
new_mu_nodes = mu_nodes;
new_lambda_edges = lambda_edges;
new_mu_edges = mu_edges;


dual = [];
primal = [];

x_nodes = x_nodes + 1/d;


% i = 0;
for i= 1:1000
    fprintf('\n%d.\n',i);
    %tStart = tic; 
    
    %% Update x_nodes
    
    lambda_nodes = new_lambda_nodes;
    mu_nodes = new_mu_nodes;
    lambda_edges = new_lambda_edges;
    mu_edges = new_mu_edges;
    
    %tic
    x_nodes = zeros(d,nNodes);
    [~,ind] = min(theta_nodes + lambda_nodes + mu_nodes,[],1);
    idx = sub2ind(size(x_nodes),ind,1:nNodes);
    x_nodes(idx) = 1;
    %toc
    
    %% Solve the edge slaves
    %tic
    
    
    
    
    D_lambda_edges = lambda_edges(kron((1:d)',ones(d,1)),:);
    C_mu_edges = repmat(mu_edges,d,1); 
    
%     x_edges = (theta_edges + D_lambda_edges + C_mu_edges) < 0;  
    
    x_edges = zeros(d*d,nEdges);    
    [~,ind2] = min(theta_edges + D_lambda_edges + C_mu_edges,[],1);
    idx2 = sub2ind(size(x_edges),ind2,1:nEdges);
    x_edges(idx2) = 1;
    
    %% Compute necessary quantities
    
    %%% Compute y_pq = D*x_pq and z_pq = C*x_pq
    M =  reshape(x_edges,[d,d,nEdges]);
    y_edges = squeeze(sum(M,1));
    z_edges = squeeze(sum(M,2));   
    
    
    x_nodes_vec = reshape(x_nodes,d*nNodes,1);
    y_edges_vec = reshape(y_edges,d*nEdges,1);
    z_edges_vec = reshape(z_edges,d*nEdges,1);
    
              
    dual_current = x_nodes_vec'*reshape(theta_nodes + lambda_nodes + mu_nodes,d*nNodes,1) + reshape(x_edges,d*d*nEdges,1)'*reshape(theta_edges + D_lambda_edges + C_mu_edges,d*d*nEdges,1);
    
    primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    
    dual(i) = dual_current;
    primal(i) = primal_current;
    
    %norm_of_gradient = 2*norm(x_nodes_vec)^2 + norm(y_edges_vec)^2 + norm(z_edges_vec)^2;     
      
    alpha = (dual_best - dual_current)/norm_of_gradient;  
    
    %primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    fprintf('Dual: %d, Primal: %d',dual_current,primal_current);
    
    if(dual_best - dual_current < 1)
        break;
    end
    
    
    %% Update the parameters    
    
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
    
    %theta_nodes = theta_nodes + alpha*(2*x_nodes - avg_out - avg_in);
    new_lambda_nodes    = lambda_nodes + alpha*(x_nodes - avg_out);
    new_mu_nodes        = mu_nodes + alpha*(x_nodes - avg_in);
    
    %toc
    
    for e = 1:nEdges
        p = edgeStruct.edgeEnds(e,1);
        q = edgeStruct.edgeEnds(e,2);
        delta_out(:,e) = y_edges(:,e) - avg_out(:,p);
        delta_in(:,e) = z_edges(:,e) - avg_in(:,q);
    end
    
    %%% Update theta_pq = theta_pq + alpha*D'*delta_out + alpha*C'*delta_in
    %theta_edges = theta_edges + alpha*delta_out(kron((1:d)',ones(d,1)),:) + alpha*repmat(delta_in,d,1);
    new_lambda_edges    = lambda_edges + alpha*delta_out;
    new_mu_edges        = mu_edges + alpha*delta_in;
    
    
    
    %% Check the projection properties
%     for p = 1:size(x_nodes,2)
%         edges1 = edges_out{p};
%         if norm(lambda_nodes(:,p) + sum(lambda_edges(:,edges1),2)) > 0.000005
%             fprintf('Lambda. p = %d\n',p);
%         end
% 
%         edges2 = edges_in{p};
%         if norm(mu_nodes(:,p) + sum(mu_edges(:,edges2),2)) > 0.000005
%             fprintf('Mu. p = %d\n',p);
%         end
%     end
    
    %alpha = alpha*2;
    
    
% 
% %     primal_current = maxVal*computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
% %     fprintf('Primal: %d, ',primal_current);
end


Energy_DD = -UGM_LogConfigurationPotential(ind,nodePot2,edgePot2,edgeStruct.edgeEnds);

figure;imagesc(reshape(ind,nRows,nCols)');title(sprintf('DD (Energy = %f)',Energy_DD));colormap gray;pause(1);
print('-depsc','results');


% save('dual','dual');
% save('primal','primal');

i_primal = 1:length(primal);
i_dual = 1:length(dual);


figure;plot(i_primal, primal,i_dual, dual);
xlim([0,4201]);
xlabel('Iterations')
ylabel('Energy');
hleg1 = legend('Primal','Dual');
title('Dual decomposition');
print('-depsc','energy');
