## Empirical Bayes Approach
Obj.eb <- MakeADFun(data=Data.tmb, parameters=Params.tmb,
                    DLL="U_EB", random=c("S_j"), silent=T)
Opt.eb <- optim(par=Obj.eb$par, fn=Obj.eb$fn, gr=Obj.eb$gr, method="BFGS",
                hessian=T)

#Opt.eb$par  # 0.688, 2.30, -3.24
cov.eb <- solve(Opt.eb$hessian)