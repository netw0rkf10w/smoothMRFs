%% Load data

% clear all
% close all
% 
% [nodePot2,edgePot2,edgeStruct,nRows,nCols,maxVal] = middlebury_loadDataUGM('tsukuba');
% [nNodes,nStates] = size(nodePot2);
% 
% 
% nodePot = -log(nodePot2);
% edgePot = -log(edgePot2);
% 
% d = nStates;
% nNeighbors = 4;
% one = ones(1,d);

theta_nodes = nodePot';
theta_edges = zeros(d*d,nNeighbors,nNodes);
delta_edges = zeros(d,nNeighbors,nNodes);
x_nodes = zeros(d,nNodes);

y_average = zeros(size(theta_nodes));

%% 
for p = 1:nNodes
    for q = 1:nNeighbors
        theta_edges(:,q,p) = 0.5*reshape(edgePot(:,:,15),d*d,1);
    end
end

%% Decode by DD

dual_current = 0;
dual_bestsofar = 0;
primal_current = 0;

alpha = 1;

%while(1)
for i=1:200
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
    x_edges = theta_edges < 0;    
    
    
%     for i=1:nNodes
%         for j=1:nNeighbors
%              y_edges(:,j,i) = (onevec*reshape(x_edges(:,j,i),d,d))';
%         end
%     end
    
    y_edges = squeeze(sum(reshape(x_edges,[d,d,nNeighbors,nNodes])));
    
    %% Update the parameters    
    y_average = squeeze(mean(y_edges,2));
    
    y_average = (nNeighbors*y_average + x_nodes)/(nNeighbors + 1);
    
       
    theta_nodes = theta_nodes + alpha*(x_nodes - y_average);
    
    %toc
    
%     tic    
    delta_edges = alpha*bsxfun(@minus,y_edges,permute(y_average,[1 3 2]));
    theta_edges = theta_edges + delta_edges(kron((1:d)',ones(d,1)),:,:);
%     toc
    
%     tic
%     for p = 1:nNodes   
%         for q=1:nNeighbors
% %             delta_edges(:,q,p) = alpha*(y_edges(:,q,p) - y_average(:,p)); 
% %             theta_edges(:,q,p) = theta_edges(:,q,p) + reshape(repmat(delta_edges(:,q,p)',d,1),[],1);
%             delta = alpha*(y_edges(:,q,p) - y_average(:,p)); 
%             theta_edges(:,q,p) = theta_edges(:,q,p) + reshape(repmat(delta',d,1),[],1);
% 
%         end
%     end
%     toc
%     tic
%     energy = reshape(x_nodes,d*nNodes,1)'*reshape(theta_nodes,d*nNodes,1) + reshape(x_edges,d*d*nNeighbors*nNodes,1)'*reshape(theta_edges,d*d*nNeighbors*nNodes,1);
%     energy 
%     toc
    %tElapsed = toc(tStart)
    
    x_nodes_vec = reshape(x_nodes,d*nNodes,1);
    y_edges_vec = reshape(y_edges,d*nNeighbors*nNodes,1);
    norm_of_gradient = norm(x_nodes_vec)^2 + norm(y_edges_vec)^2;    
       
    dual_current = x_nodes_vec'*reshape(theta_nodes,d*nNodes,1) + reshape(x_edges,d*d*nNeighbors*nNodes,1)'*reshape(theta_edges,d*d*nNeighbors*nNodes,1);
    dual_current = maxVal*dual_current;
    if dual_bestsofar < dual_current
        dual_bestsofar = dual_current;
    end
    
    alpha = (dual_bestsofar + 100 - dual_current)/norm_of_gradient;
    
    primal_current = maxVal*computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    fprintf('Dual: %d, Best dual: %d, Primal: %d, ',dual_current,dual_bestsofar,primal_current);
end


Energy_DD = -maxVal*UGM_LogConfigurationPotential(ind,nodePot2,edgePot2,edgeStruct.edgeEnds);

figure;imagesc(reshape(ind,nRows,nCols)');title(sprintf('DD (Energy = %f)',Energy_DD));colormap gray;pause(1);



% x_edges = randn(d*d,nNeighbors,nNodes);
% 
% y0 = zeros(d,nNeighbors,nNodes);
% y1 = zeros(d,nNeighbors,nNodes);
% 
% % tic
% % for i=1:nNodes
% %     for j=1:nNeighbors
% %          y0(:,j,i) = (one*reshape(x_edges(:,j,i),d,d))';
% %     end
% % end
% % toc
% 
% % tic
% % for i=1:nNodes
% %     for j=1:nNeighbors
% %          y1(:,j,i) = sum(reshape(x_edges(:,j,i),d,d))';
% %     end
% % end
% % toc
% 
% tic
% %M=reshape(M,d,d,m,n);
% y4= squeeze(sum(reshape(x_edges,[d,d,nNeighbors,nNodes])));
% toc
% 
% tic
% y2 = reshape(sum(reshape(x_edges,d,[])),d,nNeighbors,nNodes);
% toc
% 
% tic
% y3 = squeeze(sum(reshape(x_edges,[d,d,size(x_edges,2),size(x_edges,3)])));
% toc
% 




% 
% 
% 
% %% Decode by ignoring edges
% [junk,yInd] = max(nodePot,[],2);
% Energy_Ind = -maxVal*UGM_LogConfigurationPotential(yInd,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yInd,nRows,nCols)');title(sprintf('Independent (Energy = %f)',Energy_Ind));colormap gray;pause(1);
% 
% %% Decode with ICM
% yICM = UGM_Decode_ICM(nodePot,edgePot,edgeStruct);
% Energy_ICM = -maxVal*UGM_LogConfigurationPotential(yICM,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yICM,nRows,nCols)');title(sprintf('ICM (Energy = %f)',Energy_ICM));colormap gray;pause(1);
% 
% %% Decode with Alpha-Beta
% ySwap = UGM_Decode_AlphaBetaSwap(nodePot,edgePot,edgeStruct,@UGM_Decode_GraphCut,ones(nNodes,1));
% Energy_Swap = -maxVal*UGM_LogConfigurationPotential(ySwap,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(ySwap,nRows,nCols)');title(sprintf('Alpha-Beta Swap (Energy = %f)',Energy_Swap));colormap gray;pause(1);
% 
% %% Decode with Alpha-Expansion
% yExpand = UGM_Decode_AlphaExpansion(nodePot,edgePot,edgeStruct,@UGM_Decode_GraphCut,ones(nNodes,1));
% Energy_Expand = -maxVal*UGM_LogConfigurationPotential(yExpand,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yExpand,nRows,nCols)');title(sprintf('Alpha-Expansion (Energy = %f)',Energy_Expand));colormap gray;pause(1);
% 
% %% Decode with Alpha-Expansion Beta-Shrink (beta = alpha+1)
% betaSelect = 5;
% yAlphaBeta = UGM_Decode_AlphaExpansionBetaShrink(nodePot,edgePot,edgeStruct,@UGM_Decode_GraphCut,betaSelect,ones(nNodes,1));
% Energy_AlphaBeta = -maxVal*UGM_LogConfigurationPotential(yAlphaBeta,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yAlphaBeta,nRows,nCols)');title(sprintf('Alpha-Expanion Beta-Shrink with beta=alpha+1 (Energy = %f)',Energy_AlphaBeta));colormap gray;pause(1);
% 
% %% Decode with Alpha-Expansion Beta-Shrink
% betaSelect = 3;
% yAlphaBeta = UGM_Decode_AlphaExpansionBetaShrink(nodePot,edgePot,edgeStruct,@UGM_Decode_GraphCut,betaSelect,ones(nNodes,1));
% Energy_AlphaBeta = -maxVal*UGM_LogConfigurationPotential(yAlphaBeta,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yAlphaBeta,nRows,nCols)');title(sprintf('Alpha-Expanion Beta-Shrink with all beta (Energy = %f)',Energy_AlphaBeta));colormap gray



% d = 100;
% onevec = ones(1,d);
% D = kron(eye(d,d),onevec);
% Dsparse = sparse(D);
% 
% x = randn(d^2,1);
% 
% disp('dense matrix multiply')
% tic
% aa = D*x;
% toc
% 
% disp('sparse matrix multiply')
% tic
% bb = Dsparse*x;
% toc
% 
% disp('for loop version')
% tic
% cc = zeros(d,1);
% ind = 1;
% for kk=1:d
%     for jj=1:d
%         cc(kk) = cc(kk) + x(ind);
%         ind = ind + 1;
%     end
% end
% toc
% 
% disp('tensorized version')
% tic
% dd = sum(reshape(x,d,d))';
% toc
% 
% disp('ee version')
% tic
% ee = (onevec*reshape(x,d,d))';
% toc
% 
% if (norm(aa - bb) > 1e-9 || norm(bb - cc) > 1e-9 ...
%         || norm(cc - dd) > 1e-9)
%     disp('error: different methods give different results')
% end
% 
% norm(ee - dd)