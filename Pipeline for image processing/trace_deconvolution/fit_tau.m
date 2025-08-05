% 定义函数
% 定义函数
% 定义函数
func = @(param,x) max(0, param(1)*(exp(log(exp(-1./param(2)))*(x-param(4))) - exp(log(exp(-1./param(3)))*(x-param(4)))) / (exp(-1./param(2)) -exp(-1./param(3))));


ydata = double(C_Raw_g_flt(5052,2627:2709)); 
xdata = 1:length(ydata); 
plot(xdata,ydata)
A_guess=max(ydata);
% 初始参数猜测
initial_guess = [A_guess,15,0.5,5]; % 这里你需要根据你的问题来猜测一个初始值



% 使用 lsqcurvefit 函数来找到最佳参数
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
[param,resnorm,residual,exitflag,output] = lsqcurvefit(func,initial_guess,xdata,ydata,[],[],options);

% 输出结果
A = param(1)
taus1 = param(2)
taus2 = param(3)
hold on
plot(xdata,func(param,xdata),'r')
