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
% 
% 
% nodePot = -log(nodePot2);
% edgePot = -log(edgePot2);


d = nStates;
n_Ch = 2;
n_Pa = 2;

one = ones(1,d);

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

%% Decode by DD

dual_current    = 0;
dual_bestsofar  = 0;
primal_current  = 0;

alpha = 0.1;

dual_best = 370803;

for i=1:2
    fprintf('\n%d.\n',i);
    %tStart = tic; 
    %% Solve the node slaves
    
    %tic
    x_nodes = zeros(d,nNodes);
    [~,ind] = min(theta_nodes,[],1);
    idx = sub2ind(size(x_nodes),ind,1:nNodes);
    x_nodes(idx) = 1;
    %toc
    
    %% Solve the edge slaves
    %tic
    %x_edges = theta_edges < 0;  
    
    x_edges = zeros(d*d,nEdges);
    [~,ind2] = min(theta_edges,[],1);
    idx2 = sub2ind(size(x_edges),ind2,1:nEdges);
    x_edges(idx2) = 1;
    
    %% Compute necessary quantities
    
    %%% Compute y_pq = D*x_pq and z_pq = C*x_pq
    M =  reshape(x_edges,[d,d,nEdges]);
    y_edges = squeeze(sum(M,1));
    z_edges = squeeze(sum(M,2));
    
    %%% Compute the averages
    get_avg(avg_out,avg_in,x_nodes,y_edges,z_edges,edges_out,edges_in);
    
    
    %% Update the parameters    
    
    theta_nodes = theta_nodes + alpha*(2*x_nodes - avg_out - avg_in);
    
    %toc
    
    for e = 1:nEdges
        p = edgeStruct.edgeEnds(e,1);
        q = edgeStruct.edgeEnds(e,2);
        delta_out(:,e) = y_edges(:,e) - avg_out(:,p);
        delta_in(:,e) = z_edges(:,e) - avg_in(:,q);
    end
    
    %%% Update theta_pq = theta_pq + alpha*D'*delta_out + alpha*C'*delta_in
    theta_edges = theta_edges + alpha*delta_out(kron((1:d)',ones(d,1)),:) + alpha*repmat(delta_in,d,1);
    
    x_nodes_vec = reshape(x_nodes,d*nNodes,1);
    y_edges_vec = reshape(y_edges,d*nEdges,1);
    z_edges_vec = reshape(z_edges,d*nEdges,1);
    
              
    dual_current = x_nodes_vec'*reshape(theta_nodes,d*nNodes,1) + reshape(x_edges,d*d*nEdges,1)'*reshape(theta_edges,d*d*nEdges,1);
    
    dual_current = maxVal*dual_current;
    if dual_bestsofar < dual_current
        dual_bestsofar = dual_current;
    end
    
    norm_of_gradient = 2*norm(x_nodes_vec)^2 + norm(y_edges_vec)^2 + norm(z_edges_vec)^2; 
    alpha = (dual_bestsofar + 100 - dual_current)/norm_of_gradient;
    
    
    primal_current = maxVal*computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    fprintf('Dual: %d, Best dual: %d, Primal: %d, ',dual_current,dual_bestsofar,primal_current);

%     primal_current = maxVal*computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
%     fprintf('Primal: %d, ',primal_current);
end


Energy_DD = -maxVal*UGM_LogConfigurationPotential(ind,nodePot2,edgePot2,edgeStruct.edgeEnds);

figure;imagesc(reshape(ind,nRows,nCols)');title(sprintf('DD (Energy = %f)',Energy_DD));colormap gray;pause(1);

