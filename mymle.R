library('nloptr')



mymle<-function(logdensity,x,del,param0,args){

  objfun <- function(param){
    # Set the objective function
    objfun <- -logdensity2loglik(logdensity,x,del,param,args)
  }  

# Optimize
# nloptr.print.options()  


res<- nloptr( x0=param0, 
              eval_f=objfun, 
              eval_g_ineq=eval_g_ineq,
              lb = args$l, 
              ub = args$u, 
              opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = args$maxiter, "xtol_rel" = args$eps, "print_level"=args$print_level))

return(res)
# fprintf('\nComputing Variance...\n')
# % Compute the var-cov matrix
# % The derivatives are computed numerically. The Hessian may not be
# % positive definite. We report the inverse[Infomation], as well as the
# % robust sandwich matrix.
# H = hessian(objfun,param_mle);
# InfoMatrix = logdensity2info(logdensity,x,del,param_mle);
# Variance = inv(InfoMatrix);
# invH = inv(H);
# Variance_Robust = invH * InfoMatrix * invH';
# fprintf('Estimation finished.\n\n')
# 
# % Output
# output.param = param_mle;
# output.variance = Variance;
# output.variance_robust = Variance_Robust;
# output.se = sqrt(diag(Variance));
# output.se_robust = sqrt(diag(Variance_Robust));
# output.objfun = objfun;
# output.loglik = -fval;
# output.exitflag = exitflag;
# 
# % Printing the MLE result
# disp(repmat('*',1,50))
# fprintf('%10s%10s%10s\n','COEF','SE','Robust SE');
# for k=1:length(param_mle)
# fprintf('%10.4f%10.4f%10.4f\n',output.param(k),output.se(k),output.se_robust(k));
# end
# disp(repmat('*',1,50))
# disp(output.exitmsg);
# 
# % remove the search path
# rmpath(genpath(cd));
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Helper functions
# 
# function output = logdensity2info(logdensity,x,del,param)
# % Compute the information matrix = Sum_{i=1}^n Score_i' Score_i
# n = length(x) - 1;
# output = 0;
# for i=1:n
#     tmpfun = @(theta)(logdensity(x(i+1,:),x(i,:),del,theta));
#     score_i = gradest(tmpfun,param); % This is a row vector
#     output = output + score_i' * score_i;
#     end
}