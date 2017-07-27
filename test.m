
% d = 20;
% m = 30;
% n = 40;
% 
% A = randn(d,m,n);

% d = 2; n = 4;
% A = zeros(d,d,n);
% A(:,:,1) = [1 5; 2 7];
% A(:,:,2) = [-2 3; 0 5];
% A(:,:,3) = [1 -1; 2 1];
% A(:,:,4) = [6 8; 1 7];

% B = [1 2 3 4; -1 -2 -3 -4];

% tic
% B = A(kron((1:d)',ones(d,1)),:,:);
% toc
% 
% tic
% C = reshape(permute(repmat(permute(A,[1 4 2 3]),[1 d]),[ 2 1 3 4]),d*d,m,n);
% toc
% 
% 
% norm(reshape(B,d*d*m,n)-reshape(C,d*d*m,n))

% tic
% C1 = zeros(size(A));
% for  i = 1:n
%     for j = 1:m
%         C1(:,j,i) = A(:,j,i) - B(:,i);
%     end
% end
% toc



% tic
% C2 = bsxfun(@minus,A,permute(B,[1 3 2]));
% toc
% 
% tic
% B1 = repmat(reshape(B,[d 1 n]),[1 m 1]);
% C3 = A - B1;
% toc
% 
% norm(reshape(C3,d*m,n)-reshape(C2,d*m,n))


% X = [10 ; 20 ; 30 ; 40];
% d = 3 ;
% Y = X(ceil((1:d*numel(X))/d));

% d = 1000;
% X = randn(d,1);
% 
% tic
% Y = X(ceil((1:d*numel(X))/d));
% toc
% 
% tic
% Z = reshape(repmat(X',d,1),[],1);
% toc


% d= 3;
% N = 5;
% 
% ind = [2 1 3 3 1];
% X = zeros(d,N);
% 
% idx = sub2ind(size(X),ind,1:N);
% X(idx) = 1;
% X


% d = 3; n = 4;
% A = randn(d,n);
% repmat(A,d,1)
% 
% A(kron((1:d)',ones(d,1)),:)
% reshape(A(ceil(1/d:1/d:numel(A))), d^2,n)


% betaSelect = 3;
% yAlphaBeta = UGM_Decode_AlphaExpansionBetaShrink(nodePot2,edgePot2,edgeStruct,@UGM_Decode_GraphCut,betaSelect,ones(nNodes,1));
% Energy_AlphaBeta = -maxVal*UGM_LogConfigurationPotential(yAlphaBeta,nodePot2,edgePot2,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yAlphaBeta,nRows,nCols)');title(sprintf('Alpha-Expanion Beta-Shrink with all beta (Energy = %f)',Energy_AlphaBeta));colormap gray


% 
% 
% %% Check the projection properties
% for p = 1:size(x_nodes,2)
%     edges1 = edges_out{p};
%     if norm(lambda_nodes(:,p) + sum(lambda_edges(:,edges1),2)) > 0.000005
%         fprintf('Lambda. p = %d\n',p);
%     end
% 
%     edges2 = edges_in{p};
%     if norm(mu_nodes(:,p) + sum(mu_edges(:,edges2),2)) > 0.000005
%         fprintf('Mu. p = %d\n',p);
%     end
% end
    

% tt= [];
% i=0;
% 
% for x=0.1:0.01:20
%     for y=0.1:0.01:20
%         i=i+1;
%         tt(i) = (x+y)^2+(x-1)^2/y^2+(y-1)^2/x^2;
%     end
% end
% [y,I] = min(tt);


% e = ones(3,1);
% 
% D = blkdiag(e',e',e');
% C = repmat(diag(e),1,3);
% 
% A = D'*D + C'*C;

% d = 3; n = 5;
% 
% one = ones(d,1);
% A = kron(eye(d),one');
% B = repmat(diag(one),1,d);
% M = A'*A + B'*B;
% 
% % tic
% % S = d^2*n;
% % Q1 = sparse(S, S); %'// preallocates with zeros
% % for b = 0:d^2:S-1
% %     Q1(b+(1:d^2),b+(1:d^2)) = M;
% % end
% % toc
% % 
% % figure,
% % spy(Q1)
% 
% tic
% MSparse = sparse(M);
% % replicate: please don't test with n=220512 directly
% Mblk = repmat({MSparse}, 1, n);
% % convert do block diagonal
% Q = blkdiag(Mblk{:});
% toc
% 
% % %figure,
% % %spy(Q)
% % 
% % 
% % %get the directory of your input files:
% % pathname = fileparts('/home/khue/');
% % %use that when you save
% % matfile = fullfile(pathname, 'matQ.mat');
% % 
% % save(matfile, 'Q', '-v7.3');
% % 
% % %load(matfile);


% load('1_primal');
% load('1_dual');
% 
% primal_subgrad = primal;
% dual_subgrad = dual; 
% 
% load('1_primal_smooth1');
% load('1_dual_smooth1');
% primal_smooth1 = primal - 0.01*(nNodes+2*nEdges)*log(d);
% dual_smooth1 = dual-0.01*(nNodes+2*nEdges)*log(d); 
% 
% load('1_primal_smooth3');
% % load('1_dual_smooth1');
% primal_smooth3 = primal;
% 
% 
% % primal_smooth = primal-mu*(nNodes+2*nEdges)*log(d);
% % dual_smooth = dual-mu*(nNodes+2*nEdges)*log(d); 
% 
% % save('0_primal_subgrad','primal_subgrad');
% % save('0_dual_subgrad','dual_subgrad');
% 
% 
% rsub = 1:150;
% 
% l = min([length(primal_subgrad) length(primal_smooth1) length(primal_smooth3)]);
% 
% rr=1:l;
% 
% dual_optimal = 484.7098 + zeros(1,length(rr));
% 
% 
% figure;plot(rr,dual_optimal,'r',rr, primal_subgrad(rr),'b',rr, dual_subgrad(rr),'b--',rr, primal_smooth1(rr),'g',rr,primal_smooth3(rr),'r');
% %figure;plot(rr, primal_subgrad(rr),'b',rr, dual_subgrad(rr),'b--',rr, primal_smooth(rr),'r',rr, dual_smooth(rr),'r--');
% xlim([0,length(rr)]);
% xlabel('Iterations')
% ylabel('Energy');
% hleg1 = legend('Optimal','Subgradient','Dual subgradient','Smoothed, mu=0.01','Smoothed, mu=0.001');
% title('Dual decomposition');
% print('-depsc','energy');


% tt = [-1 1.2 0.5 -0.5; 0.1 -0.1 15 -2]
% max(tt,0)


% tt = randn(d*d,nEdges);
% ss = randn(d*d,nEdges);
% 
% tic
% s1 = reshape(tt,d*d*nEdges,1)'*reshape(ss,d*d*nEdges,1);
% 
% 
% 
% 
% toc
% 
% tic
% s2 = sum(sum(tt.*ss));
% toc


%% Decode with Alpha-Expansion
% yExpand = UGM_Decode_AlphaExpansion(nodePot2,edgePot2,edgeStruct,@UGM_Decode_GraphCut,ones(nNodes,1));
% Energy_Expand = -maxVal*UGM_LogConfigurationPotential(yExpand,nodePot,edgePot,edgeStruct.edgeEnds)
% figure;imagesc(reshape(yExpand,nRows,nCols)');title(sprintf('Alpha-Expansion (Energy = %f)',Energy_Expand));colormap gray;pause(1);



% cc = 1e4*randn(10,20);
% 
% n = size(cc,1);
% 
% dd = exp(cc);
% 
% sum(sum(isnan(bsxfun(@rdivide,dd,sum(dd)))))


%A = 1e4*randn(10,20);
c= 0.9999;
A = 1e4*[1 -0.5 -1; c 0.5 -0.5; -0.5 -0.5 -0.5]
S = sym(A,'d'); %// convert to symbolic
C = exp(S)./repmat(sum(exp(S)),size(S,1),1);

Cd = eval(C);
Cd


A2 = bsxfun(@minus,A,max(A));

B = exp(A2);
C = bsxfun(@rdivide,B,sum(B))

% x = [1 2 3;2 3 4]
% v = [2 0 1]
% z=bsxfun(@minus,x,v)

