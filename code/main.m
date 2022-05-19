%%
% This is the code for "Multi-scale hybrid spherical graphite composites：a lightweight thermal interface material with high thermal conductivity and simple processing technology"

clc;clear;close all;

Dmin = 0.58; % Dmin;
Dmax = 70;   % Dmax;

%% 0. Read the data (读取数据)
% The particle size distribution data of graphite were recorded in data.xlsx 
% For each data: 
% Column one is the Particle size; 第一列是粒径
% Column two is the Cumulative volume fraction; 第二列是累积体积分数
% Column three is the Volume fraction; 第三列是体积分数
data0 = xlsread('data.xlsx');
Particle_size = data0(:,1);
[~,n] = size(data0);   
Cumulative_volume_fraction = data0(:,2:3:n-1);  % Extract cumulative volume fraction of powders
Volume_fraction = data0(:,3:3:n); % Extract volume fraction of powders
N = n/3; % There are total N powder

%% 1. Data processing: Solve the fitting problem of lognormal distribution curve (数据处理: 解决对数正态分布曲线的拟合问题)

fprintf('1. Data processing : the parameters of volume fraction are fitted by nonlinear least squares');
fprintf('\n（利用非线性最小二乘拟合体积分数的参数）');

for i = 1:N
    fprintf('\n -----------------------------------------------');
    fprintf('\n Fitting the %2d-th curve （曲线拟合第 %2d 条） \n',i, i);
    data_temp = Cumulative_volume_fraction(:,i);
    index = ( data_temp ~= 0 & data_temp ~= 100 );
    data2 = Volume_fraction(:,i);
    x = [Particle_size(find(index ~= 0,1) - 1);Particle_size(index)];
    y = 0.01*data_temp(index);
    z = 0.01*data2(index);
    [mu,sigma,~] = LM_1(x,y);          % fitting
    
    f = @(x) exp(-(log(x)-mu).^2/(2*sigma.^2))./(sqrt(2*pi)*sigma.*x);
    fprintf('\n mu = %4f',mu);          % The parameter mu
    fprintf('\n sigma = %4f',sigma);    % The parameter sigma
    fprintf('\n Expectation = %4f',exp(mu+sigma^2/2));
    fprintf('\n Std = %4f',sqrt(exp(2*mu+sigma^2)*(exp(sigma^2)-1)));
    Mu(i) = mu;
    Sigma(i) = sigma;
end


%% 2. Mathematical modeling (数学建模)
% solving parameters H and f（求解参数 H 和 f）
A = zeros(N,N);
q = zeros(N,1);
n = 0.37;
U = @(t) (t.^n - Dmin^n)./(Dmax^n - Dmin^n);   % U is Dinger-Funk eqaution
U2 = @(t) U(t).^2;
C = integral(U2,Dmin,Dmax);
fprintf('\n -----------------------------------------------');

fprintf('\n\n 2. Mathematical modeling（数学建模）:\n Calculating the symmetric matrix H and vector f ...');
%fprintf('\n 计算对称矩阵H和向量f...');
for i = 1:N
    mu = Mu(i);
    sigma = Sigma(i);
    f = @(x) exp(-(log(x)-mu).^2./(2*sigma.^2))./(x.*sqrt(2*pi)*sigma);
    F = @(y) Simpson(f,Dmin,y,1.0e-3).^2;                             % F is calculated by Simpson Formulation
    A(i,i) = compositeTrapezoidintegrationrule(F,Dmin,Dmax,1.0e-2);   % A is calculated by composite Trapezoid Formulation
    fprintf('\n H(%d,%d) = %4f',i,i,2*A(i,i));
    for j = (i+1):N
        mu2 = Mu(j);
        sigma2 = Sigma(j);
        if(i ~= j)
            f2 = @(x) exp(-(log(x)-mu2).^2./(2*sigma2.^2))./(x.*sqrt(2*pi)*sigma2);
            F2 = @(y) Simpson(f2,Dmin,y,1.0e-3).*Simpson(f,Dmin,y,1.0e-3);
            A(i,j) = compositeTrapezoidintegrationrule(F2,Dmin,Dmax,1.0e-2);
            A(j,i) = A(i,j);
            fprintf('\n H(%d,%d) = %4f',i,j,2*A(i,j));
        end
    end
    F3 = @(y) Simpson(f,Dmin,y,1.0e-3).*U(y);
    q(i) = compositeTrapezoidintegrationrule(F3,Dmin,Dmax,1.0e-2);
    fprintf('\n f(%d) = %4f',i,-2*q(i));
end

%% 3. Model Solving (模型求解)
% Calling Quadratic programming API "quadprog()" （调用二次规划程序）
if(min(eig(A)) > 0)
    H = 2*A;
    f = -2*q;
    % min 1/2 *t'*H*t + f'*t
    [t,fval,exitflag] = quadprog(H,f,[],[],ones(1,N),1,zeros(N,1));
end
%}
fprintf("\n\n3. Model Solving（模型求解）:\nThe best formulation is:");
fprintf("\n Mu           Sigma      proportion(percent)");
fprintf('\n -----------------------------------------------');
for i = 1:N
    fprintf('\n %4f     %4f     %4.2f %%',Mu(i), Sigma(i), 100 * t(i));
end


