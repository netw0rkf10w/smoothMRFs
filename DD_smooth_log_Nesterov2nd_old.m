
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


% %% The smoothing constant
% mu = 0.01;
% 
% %% Compute the Lipschitz constant L = norm(B)^2/mu
% %L = 1/mu*2*nEdges*d*(d+1);
% L = 50;
% 
% 
% theta_nodes = nodePot';
% theta_edges = reshape(edgePot,d*d,nEdges);
% 
% x_nodes = zeros(d,nNodes);
% x_edges = zeros(d*d,nEdges);
% 
% y_edges = zeros(d,nEdges);
% z_edges = zeros(d,nEdges);
% 
% avg_out = zeros(d,nNodes);
% avg_in = zeros(d,nNodes);
% 
% delta_out = zeros(d,nEdges);
% delta_in = zeros(d,nEdges);
% 
% 
% lambda_edges = zeros(d,nEdges);
% nu_edges = zeros(d,nEdges);
% t_lambda = zeros(d,nEdges);
% t_nu = zeros(d,nEdges);
% 
% g_lambda = zeros(d,nEdges);
% g_nu = zeros(d,nEdges);
% 
% 
% %% Decode by DD
% 
% dual_current    = 0;
% dual_bestsofar  = 0;
% primal_current  = 0;
% 
% dual_best = 484.7098;
% 
% 
% 
% 
% dual = [];
% primal = [];
% 
% alpha = 1;


% i = 0;
for i= 102:151
    fprintf('\n%d.\n',i);
    %tStart = tic; 
    %% Solve the node slaves 
    
    
    %% Compute the gradient by solving the node and edge subproblems
    %fprintf('Solve node slaves... \n')
    for p = 1:nNodes
        %%% Compute y_p and z_p
        edges1 = edges_out{p};
        edges2 = edges_in{p};
        lambda_p = -sum(lambda_edges(:,edges1),2);        
        nu_p = -sum(nu_edges(:,edges2),2);
        
        a = 1/mu*(theta_nodes(:,p) + lambda_p + nu_p) -2/d;
        k = (sum(a)+2)/d;
        
        x_nodes(:,p) = exp(-a)/sum(exp(-a));
    end   
    
    %% Solve the edge slaves
    %fprintf('Solve edge slaves... \n')
    %%% Compute D'*lambda_pq + C'*nu_pq
    D_lambda_edges = lambda_edges(kron((1:d)',ones(d,1)),:);
    C_nu_edges = repmat(nu_edges,d,1); 
    
        
    b = 1/mu*(theta_edges + D_lambda_edges + C_nu_edges);
    x_edges = exp(-b); 
    %%% Divide each column by its sum
    x_edges = bsxfun(@rdivide,x_edges,sum(x_edges));    
    
    
    
    %% Compute necessary quantities
    
    %%% Compute y_pq = D*x_pq and z_pq = C*x_pq
    M =  reshape(x_edges,[d,d,nEdges]);
    y_edges = squeeze(sum(M,1));
    z_edges = squeeze(sum(M,2));   
    
    %% Compute the gradient g(lambda,nu)
    for e = 1:nEdges
        p = edgeStruct.edgeEnds(e,1);
        q = edgeStruct.edgeEnds(e,2);
        g_lambda(:,e) = y_edges(:,e) - x_nodes(:,p);
        g_nu(:,e) = z_edges(:,e) - x_nodes(:,q);
    end    
    
    
           
    dual_current = reshape(theta_nodes,d*nNodes,1)'*reshape(x_nodes,d*nNodes,1)...
        + reshape(theta_edges,d*d*nEdges,1)'*reshape(x_edges,d*d*nEdges,1)...
        + reshape(lambda_edges,d*nEdges,1)'*reshape(g_lambda,d*nEdges,1)... 
        + reshape(nu_edges,d*nEdges,1)'*reshape(g_nu,d*nEdges,1);
    
    [~,ind] = max(x_nodes,[],1);
    primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    
    dual(i) = dual_current;
    primal(i) = primal_current;
    
    
    fprintf('Dual: %d, Primal: %d',dual_current,primal_current);
    
%     if(dual_best - dual_current < 1)
%         break;
%     end
    
    
    %% Update the parameters    
    
    t_lambda_new = lambda_edges + 1/L*g_lambda;
    t_nu_new = nu_edges + 1/L*g_nu;
    
    
    alpha_new = 0.5*(1+sqrt(1+4*alpha^2));
    gamma = (1-alpha)/alpha_new;
    
    lambda_edges = (1-gamma)*t_lambda_new + gamma*t_lambda;
    nu_edges = (1-gamma)*t_nu_new + gamma*t_nu;    
    
    
    t_lambda = t_lambda_new;
    t_nu =t_nu_new;
    
end

[~,ind] = max(x_nodes,[],1);

Energy_DD = -UGM_LogConfigurationPotential(ind,nodePot2,edgePot2,edgeStruct.edgeEnds);

figure;imagesc(reshape(ind,nRows,nCols)');title(sprintf('DD (Energy = %f)',Energy_DD));colormap gray;pause(1);
print('-depsc','results_DD_smooth_log_Nesterov2nd_old');


% save('dual','dual');
% save('primal','primal');

figure, plot(primal);
print('-depsc','energy_DD_smooth_log_Nesterov2nd_old');

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
