function [mu,sigma,normF] = LM_1(x,y)
%  Levenberg-Marquardt
%  To solve the parameter of mu,sigma in the lognormal distribution

tol = 1.0e-7;
mu = 2; sigma = 0.6;      % Initialize mu,sigma
rho = 0.2;
sigma_1 = 0.4;
% mu_k = 1;
k = 0;
itermax = 99;
J = JF(mu,sigma,x);
Fy = F(mu,sigma,x,y);
fprintf('------------------------------------------------\n');

fprintf('iter   i      mu      sigma      norm(J^T*Fy)   1/2*norm(Fy,2)^2');
while norm(J'*Fy) >= tol && k <= itermax
    mu_k = max(0,-min(eig(J'*J))+tol);
    d = (J'*J + mu_k*eye(2))\(-J'*Fy);
    for i = 1:20
        if(1/2*norm(F(mu+rho^(i-1)*d(1),sigma+rho^(i-1)*d(2),x,y),2)^2 <= ( 1/2*norm(Fy,2)^2 + sigma_1*rho^(i-1)*((J'*Fy)'*d))  )
            break;
        end
    end
    mu = mu + rho^(i-1)*d(1);
    sigma = sigma + rho^(i-1)*d(2);
    J = JF(mu,sigma,x);
    Fy = F(mu,sigma,x,y);
    fprintf('\n%2d   %2d   %6f    %6f    %6f    %6f',k,i,mu,sigma,norm(J'*Fy),1/2*norm(Fy,2)^2);
    k = k + 1;
end
normF = 1/2*norm(Fy,2)^2;


function [Fy] = Fx(mu,sigma,a,b,y)
    % Calculate the difference between the integral of the probability density function of lognormal distribution on [a, b] and the actual cumulative volume fraction y
    % 求对数正态分布的概率密度函数在[a,b]上的积分与实际累计体积分数y的差
    f = @(t) exp(-t.^2);
    Fy = integral(f,(log(a)-mu)/(sqrt(2)*sigma),(log(b)-mu)/(sqrt(2)*sigma))/sqrt(pi) - y;
return

function [D] = DF(mu,sigma,a,b)
    % Calculate Jacobi matrix of the integral of one of the probability density functions of lognormal distribution on [a, b] with respect to the two parameters mu, sigma
    % 求其中一个对数正太分布的概率密度函数在[a,b]上的积分关于两个参数mu,sigma的Jacobi矩阵
    DF1 = (exp(-(log(a)-mu)^2/(2*sigma^2)) - exp(-(log(b)-mu)^2/(2*sigma^2)))/(sigma*sqrt(2*pi)); 
    DF2 = ((log(a)-mu)*exp(-(log(a)-mu)^2/(2*sigma^2)) - (log(b)-mu)*exp(-(log(b)-mu)^2/(2*sigma^2)))/(sqrt(2*pi)*sigma^2); 
    D = [DF1;DF2];
return

function [J] = JF(mu,sigma,x)
    % Calculate thr Jacobi matrix of F
    % 求解F的Jacobi矩阵
    J = zeros(length(x)-1,2);
    for i = 1:length(x)-1
        J(i,:) = DF(mu,sigma,x(1),x(i+1))';
    end
return

function [FY] = F(mu,sigma,x,y)
    % F = [int_{x(i)}^{x(i+1)}p_{mu,sigma}(x)dx-y_i]_{i=1}^n
    FY = zeros(length(y),1);
    for i = 1:length(y)
        FY(i) = Fx(mu,sigma,x(1),x(i+1),y(i));
    end
return