
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
% L = 1000;


%% Initialization
L_old = 100000;
beta = 1.1;
t_old = 1;

lambda_x_old = zeros(d,nEdges);
nu_x_old = zeros(d,nEdges);
lambda_y = zeros(d,nEdges);
nu_y = zeros(d,nEdges);


dual_current    = 0;
dual_bestsofar  = 0;
primal_current  = 0;

Ls = [];


for k= 1:200
    fprintf('\n%d.\n',k);
    %tStart = tic; 
    %% Solve the node slaves 
    
    L = L_old;    
    
    
    %% Compute the gradient
    
    [func_y, grad_lambda, grad_nu, ~, ~] =  compute(lambda_y, nu_y, mu, theta_nodes, theta_edges, edgeStruct, edges_out, edges_in);
     
    lambda_x = lambda_y - 1/L*grad_lambda;
    nu_x = nu_y - 1/L*grad_nu;
    
    [func_x, ~, ~, u_nodes, u_edges] =  compute(lambda_x, nu_x, mu, theta_nodes, theta_edges, edgeStruct, edges_out, edges_in);
   

    lambda_dd = lambda_x - lambda_y;
    nu_dd = nu_x - nu_y;   
    
    
    
    
    while func_x > (func_y + reshape(grad_lambda,d*nEdges,1)'*reshape(lambda_dd,d*nEdges,1)...
            + reshape(grad_nu,d*nEdges,1)'*reshape(nu_dd,d*nEdges,1)...
            + L/2*(norm(lambda_dd,'fro')^2 + norm(nu_dd,'fro')^2)),
        
        fprintf('Back-tracking\n');        
        L = beta*L;
        
        fprintf('L=%f\n',L);
        
        [func_y, grad_lambda, grad_nu, ~, ~] =  compute(lambda_y, nu_y, mu, theta_nodes, theta_edges, edgeStruct, edges_out, edges_in);
     
%         lambda_x = lambda_y - 1/L*grad_lambda;
%         nu_x = nu_y - 1/L*grad_nu;
        
        lambda_x = lambda_x_old - 1/L*grad_lambda;
        nu_x = nu_x_old - 1/L*grad_nu;

        [func_x, ~, ~, u_nodes, u_edges] =  compute(lambda_x, nu_x, mu, theta_nodes, theta_edges, edgeStruct, edges_out, edges_in);


        lambda_dd = lambda_x - lambda_y;
        nu_dd = nu_x - nu_y; 
        
    end
    
    fprintf('func_y = %f, func_x = %f\n',func_y,func_x);
    
    t = (1+ sqrt(1+4*t_old^2))/2; 
    
        
    lambda_y = lambda_x + (t_old-1)/t*(lambda_x - lambda_x_old);
    nu_y = nu_x + (t_old-1)/t*(nu_x - nu_x_old);
    
    
    %% Update
    lambda_x_old = lambda_x;
    nu_x_old = nu_x;
    t_old = t;       
    
    Ls(k) = L;
    L_old = L;
    
    [~,ind] = max(u_nodes,[],1);
    %primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds) + mu*(nNodes+2*nEdges)*log(d);
    
    primal_current = computeEnergy(ind,nodePot,edgePot,edgeStruct.edgeEnds);
    %dual(k) = dual_current;
    primal(k) = primal_current;    
    
    
    %fprintf('Dual: %d, Primal: %d',dual_current,primal_current);
    fprintf('Primal: %d',primal_current);
    
    if(primal_current - dual_best < 1)
        break;
    end
    
    
    
end



%[~,ind] = max(u_nodes,[],1);

Energv_DD = -UGM_LogConfigurationPotential(ind,nodePot2,edgePot2,edgeStruct.edgeEnds);

figure;imagesc(reshape(ind,nRows,nCols)');title(sprintf('DD (Energy = %f)',Energv_DD));colormap gray;pause(1);
print('-depsc','results_FISTA_backtracking_mu_1e5');


figure, plot(primal);
print('-depsc','energy_FISTA_backtracking_mu_1e5');

figure, plot(Ls);



%save('1_dual_FISTA_backtracking','dual');
save('1_primal_FISTA_backtracking_mu_1e5','primal');



% 
% i_primal = 1:length(primal);
% i_dual = 1:length(dual);
% 
% 
% 
% figure;plot(i_primal, primal,i_dual, dual);
% xlim([0,500]);
% xlabel('Iterations')
% ylabel('Energy');
% hleg1 = legend('Primal','Dual');
% title('Dual decomposition');
% print('-depsc','energy');

% delete(gcp('nocreate'));
% matlabpool('close');
