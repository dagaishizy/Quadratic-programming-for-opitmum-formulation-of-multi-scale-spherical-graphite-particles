This code runs in Matlab environment (本代码运行于Matlab环境)

You need to put your distribution data of powder as the form as xlsx file descriped by "0. Read the data" below, then you would get the optimum formulation by "main.m".
(你需要把粉末的分布数据排成下面叙述中 "0.Read the data" 描述的的xlsx格式，即可通过"main.m"可以得到最佳配比)

Introduction  of program "main.m" (程序"main.m"的说明)：
0. Read the data (0. 读取数据);
	The data in the xlsx file is the particle size distribution data of graphite （xlsx文件中的数据是石墨的粒度分布数据）
	Each powder has three column data （每种粉末都有三列的数据）
	1) Column one is the Particle size (um); （第一列是粒径）
	2) Column two is the Cumulative volume fraction (%); （第二列是累积体积分数）
	3) Column three is the Volume fraction (%); （第三列是体积分数）
1. Data processing: Solve the fitting problem of lognormal distribution curve (1. 数据处理: 解决对数正态分布曲线的拟合问题)
	Obtain the parameters "mu" and "sigma" of each powder, as well as the corresponding distribution "expectation", variance "std"
	（得到每种粉末的参数 mu 和 sigma ，以及对应的分布期望Expectation, 方差Std）
2. Mathematical modeling (2. 数学建模)
	Solving parameters in the quadratic programming which is H and f
	（求解二次规划参数 H 和 f）
3.  Model Solving (3. 模型求解)
	The best formulation is solved by quadratic programming. The best proportion is recorded as t. 
	(求解二次规划问题得到最佳配比t)

Remark:
	You may need to set the parameter "Dmin" and "Dmax" for your problem, In this example, we set Dmin = 0.5(um), Dmax = 70(um). 
	(你可能需要针对你的问题设置参数"Dmin"和 "Dmax"，本案例中设置Dmin = 0.5(um), Dmax = 70(um))
