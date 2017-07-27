%% Load data
% 
% % clear all
% % close all
% 
% lambda = 0.01;     
% K = 2;
% 
% 
% 
% %% Read the images
% im1 = im2double(rgb2gray(imread('tsukuba/tsukuba_left.png')));
% im2 = im2double(rgb2gray(imread('tsukuba/tsukuba_right.png')));
% % figure, imshow(im1);
% 
% %% Get sizes and dimensions 
% [h,w] = size(im1);
% 
% nStates = 16;
% dmax = nStates - 1; % dmin = 0
% 
% % Radius of the patch
% r = 15;
% 
% 
% h_crop = h - 2*r;
% w_crop = w - 2*r - dmax;
% 
% nNodes = w_crop*h_crop;
% nEdges = h_crop*(w_crop-1) + w_crop*(h_crop-1); % Grid graph
% 
% % 
% nodePot = zeros(nNodes,nStates);
% 
% 
% 
% %% Compute data terms
% 
% p=0;
% for y=(r+1):(h-r)
%     for x=(dmax+r+1):(w-r)
%         p = p+1;
%         for d = 0:dmax
%             nodePot(p,d+1) =  sum_square(y,x,y,x-d,r,im1,im2);
%             %nodePot(p,d+1) =  affine(y,x,y,x-d,r,im1,im2);
%             %nodePot(p,d+1) = abs(im1(y,x) - im2(y,x-d));
%         end
%     end    
% end
% 
% 
% 
% % save('nodePot','nodePot');
% % 
% % load('nodePot');
% 
% 
% %% Compute smoothness terms
% % edgePot = zeros(nStates,nStates,nEdges);
% 
% % [l0,l1] = meshgrid( 0:nStates-1, 0:nStates-1 );
% % dist = abs(l0-l1);
% % 
% 
% 
% edgePot = repmat(lambda*min(dist,K),[1,1,nEdges]);
% 
% %maxVal = max([nodePot(:);edgePot(:)]);
% 
% nodePot = exp(-nodePot);
% edgePot = exp(-edgePot);
% 
% 
% %% Make edgeStruct
% 
% fprintf('Making Adjacency Matrix...\n');


% % adj = zeros(nNodes,nNodes);
% % p = 0;
% % for y=1:h_crop
% %     for x=1:w_crop
% %         p = p+1;
% %         if x < w_crop
% %             adj(p,p+1) = 1;
% %             adj(p+1,p) = 1;
% %         end
% %         if y < h_crop
% %             adj(p,p+w_crop) = 1;
% %             adj(p+1,w_crop) = 1;
% %         end
% %     end
% % end
% % 
% % edgeStruct = UGM_makeEdgeStruct(adj,nStates);



% ind1full = zeros(0,1);
% ind2full = zeros(0,1);
% for c = 1:w_crop
% 	ind1 = [1:h_crop-1] + h_crop*(c-1);
% 	ind2 = [2:h_crop] + h_crop*(c-1);
% 	ind1full(end+1:end+length(ind1)) = ind1;
% 	ind2full(end+1:end+length(ind2)) = ind2;
% end
% for r = 1:h_crop
% 	ind1 = r + h_crop*([1:w_crop-1]-1);
% 	ind2 = r + h_crop*([2:w_crop]-1);
% 	ind1full(end+1:end+length(ind1)) = ind1;
% 	ind2full(end+1:end+length(ind2)) = ind2;
% end
% adj = sparse(ind1full,ind2full,1,nNodes,nNodes);
% adj = adj+adj';
% 
% fprintf('Making edgeStruct...\n');
% edgeStruct = UGM_makeEdgeStruct(adj,nStates);
% 
nRows = h_crop;
nCols = w_crop;

%% Decode by ignoring edges
[junk yInd] = max(nodePot,[],2);
figure;imagesc(reshape(yInd,nCols,nRows)');title(sprintf('Independent (Energy = %f)',Energy_Ind));colormap gray;pause(1);

% %% Decode with Alpha-Expansion
% yExpand = UGM_Decode_AlphaExpansion(nodePot,edgePot,edgeStruct,@UGM_Decode_GraphCut,ones(nNodes,1));
% %Energy_Expand = -maxVal*UGM_LogConfigurationPotential(yExpand,nodePot,edgePot,edgeStruct.edgeEnds);
% figure;imagesc(reshape(yExpand,nCols,nRows)');title(sprintf('Alpha-Expansion (Energy = %f)',Energy_Expand));colormap gray;pause(1);

%% Decode with Alpha-Expansion Beta-Shrink
% betaSelect = 3;
% yAlphaBeta = UGM_Decode_AlphaExpansionBetaShrink(nodePot,edgePot,edgeStruct,@UGM_Decode_GraphCut,betaSelect,ones(nNodes,1));
% Energy_AlphaBeta = -maxVal*UGM_LogConfigurationPotential(yAlphaBeta,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yAlphaBeta,nCols,nRows)');title(sprintf('Alpha-Expanion Beta-Shrink with all beta (Energy = %f)',Energy_AlphaBeta));colormap gray



