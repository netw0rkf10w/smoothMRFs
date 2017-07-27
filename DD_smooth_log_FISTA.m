
%% Load data

% myPool = parpool(4);

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



theta_nodes = nodePot';
theta_edges = reshape(edgePot,d*d,nEdges);

% theta_nodes = 0.0001*theta_nodes;
% theta_edges = 0.0001*theta_edges;

lambda_nodes = zeros(d,nNodes);
nu_nodes = zeros(d,nNodes);


%% Allocation
u_nodes = zeros(d,nNodes);
u_edges = zeros(d*d,nEdges);
grad_lambda_y = zeros(d,nEdges);
grad_nu_y = zeros(d,nEdges);


%% Decode by DD



dual_best = 484.7098;


dual = [];
primal = [];

stsize = [];


%% The smoothing constant
mu = 0.00001;

%% Compute the Lipschitz constant L = norm(B)^2/mu
% L = 4*(nNodes + nEdges)/mu;
L = 100000;


%% Initialization

t_old = 1;
lambda_x_old = zeros(d,nEdges);
nu_x_old = zeros(d,nEdges);
lambda_y = zeros(d,nEdges);
nu_y = zeros(d,nEdges);


dual_current    = 0;
dual_bestsofar  = 0;
primal_current  = 0;


for k= 1:201
    fprintf('\n%d.\n',k);
    
    %%
    %% Compute the gradient df(y_k)
    %%
    
    
%     parfor p = 1:nNodes
%         %%% Compute v_p and w_p
%         edges1 = edges_out{p};
%         edges2 = edges_in{p};
%         lambda_nodes(:,p) = -sum(lambda_y(:,edges1),2);        
%         nu_nodes(:,p) = -sum(nu_y(:,edges2),2);
%     end 
%     
%     u_nodes = zeros(d,nNodes);
%     [~,ind] = min(theta_nodes + lambda_nodes + nu_nodes,[],1);
%     idx = sub2ind(size(u_nodes),ind,1:nNodes);
%     u_nodes(idx) = 1;
    
    
    
    parfor p = 1:nNodes
        %%% Compute v_p and w_p
        edges1 = edges_out{p};
        edges2 = edges_in{p};
        lambda_p = -sum(lambda_y(:,edges1),2);        
        nu_p = -sum(nu_y(:,edges2),2);
        
        a = -1/mu*(theta_nodes(:,p) + lambda_p + nu_p);
        
        a = a - max(a);
        
        u_nodes(:,p) = exp(a)/sum(exp(a));
    end   
    
    if(sum(sum(isnan(u_nodes))) > 0)
        fprintf('NaN !!!!!!\n');
    end



%% Solve the edge slaves
    %fprintf('Solve edge slaves... \n')
    %%% Compute D'*lambda_pq + C'*nu_pq
    D_lambda = lambda_y(kron((1:d)',ones(d,1)),:);
    C_nu = repmat(nu_y,d,1);    
    
        

%     u_edges = zeros(d*d,nEdges);    
%     [~,ind2] = min(theta_edges + D_lambda + C_nu,[],1);
%     idx2 = sub2ind(size(u_edges),ind2,1:nEdges);
%     u_edges(idx2) = 1;

    
    b = -1/mu*(theta_edges + D_lambda + C_nu);
    b = bsxfun(@minus,b,max(b));
    
    u_edges = exp(b); 
    %%% Divide each column by its sum
    u_edges = bsxfun(@rdivide,u_edges,sum(u_edges));    

    if(sum(sum(isnan(u_edges))) > 0)
        fprintf('NaN !!!!!!\n');
    end
    
    %% Compute necessary quantities
    
    %%% Compute v_pq = D*u_pq and w_pq = C*u_pq
    M =  reshape(u_edges,[d,d,nEdges]);
    v_edges = squeeze(sum(M,1));
    w_edges = squeeze(sum(M,2));   
    
    %% Compute the gradient g(lambda,nu)
    parfor e = 1:nEdges
        p = edgeStruct.edgeEnds(e,1);
        q = edgeStruct.edgeEnds(e,2);
        grad_lambda(:,e) = -v_edges(:,e) + u_nodes(:,p);
        grad_nu(:,e) = -w_edges(:,e) + u_nodes(:,q);
    end  
    %%
    %% End of computation of the gradient
    %%
    
    
    lambda_x = lambda_y - 1/L*grad_lambda;
    nu_x = nu_y - 1/L*grad_nu;
    %% Replacing with the following, works very well but do not know why
    %% lambda_x = lambda_x_old - 1/L*grad_lambda;
    %% nu_x = nu_x_old - 1/L*grad_nu;
    
    t = (1+ sqrt(1+4*t_old^2))/2;
    
    
        
    lambda_y = lambda_x + (t_old-1)/t*(lambda_x - lambda_x_old);
    nu_y = nu_x + (t_old-1)/t*(nu_x - nu_x_old);
    
    
    %% Update
    lambda_x_old = lambda_x;
    nu_x_old = nu_x;
    t_old = t;

    
    [~,ind] = max(u_nodes,[],1);
    %primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds) + mu*(nNodes+2*nEdges)*log(d);
    
    primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    %dual(k) = dual_current;
    primal(k) = primal_current;    
    
    
    
    %fprintf('Dual: %d, Primal: %d',dual_current,primal_current);
    fprintf('Primal: %d',primal_current);
    
%     if(primal_current - dual_best < 10)
%         break;
%     end
    
    
    
end



%[~,ind] = max(u_nodes,[],1);

Energv_DD = -UGM_LogConfigurationPotential(ind,nodePot2,edgePot2,edgeStruct.edgeEnds);

figure;imagesc(reshape(ind,nRows,nCols)');title(sprintf('DD (Energy = %f)',Energv_DD));colormap gray;pause(1);
%print('-depsc','results_FISTA_15000');

% 
% 
% save('1_dual_smooth3','dual');
save('1_primal_FISTA_L15000','primal');
% 

figure,
plot(primal);
%print('-depsc','energy_FISTA_L15000');

% 
% i_primal = 1:length(primal);
% i_dual = 1:length(dual);
% 
% 
% 
% figure;plot(i_primal, primal,i_dual, dual);
% xlim([0,length(i_primal)]);
% xlabel('Iterations')
% ylabel('Energy');
% hleg1 = legend('Primal','Dual');
% title('Dual decomposition');
 



% delete(gcp('nocreate'));
