


load('1_primal');
load('1_dual');

primal_subgrad = primal;
dual_subgrad = dual; 


load('1_primal_FISTA_backtracking');
primal_FISTA_backtracking = primal;


load('1_primal_FISTA_backtracking_mu_1e5');
primal_FISTA_backtracking_mu_1e5 = primal;


% primal_smooth = primal-mu*(nNodes+2*nEdges)*log(d);
% dual_smooth = dual-mu*(nNodes+2*nEdges)*log(d); 

% save('0_primal_subgrad','primal_subgrad');
% save('0_dual_subgrad','dual_subgrad');


%rsub = 1:150;

% l = min([length(primal_subgrad) length(primal_FISTA_backtracking) length(primal_FISTA_backtracking_mu_1e5)]);

rr=1:l;



dual_optimal = 484.7098 + zeros(1,length(rr));


figure;plot(rr,dual_optimal,'k',rr, primal_subgrad(rr),'b',rr, dual_subgrad(rr),'b--',rr, primal_FISTA_backtracking(rr),'r');
%figure;plot(rr, primal_subgrad(rr),'b',rr, dual_subgrad(rr),'b--',rr, primal_smooth(rr),'r',rr, dual_smooth(rr),'r--');
xlim([0,length(rr)]);
xlabel('Iterations')
ylabel('Energy');
hleg1 = legend('Optimal','Subgradient','Dual subgradient','FISTA-backtracking, mu=0.0001','FISTA-backtracking, mu=0.00001');
title('Dual decomposition');
print('-depsc','energy_FISTA2');


