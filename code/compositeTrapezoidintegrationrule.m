function [fx]=compositeTrapezoidintegrationrule(f,a,b,h)
%   复化梯形公式
%   Input：f为要求导的函数；a,b为所求的函数的积分区间;h为步长,其中（b-a）取为h的整数倍;
%   Output：fx为f在[a,b]上的积分值
    assert(nargin==4);
    n = round((b-a)/h);
%assert((rem(n,1)==0),['更改步长h,使得(b-a)取为h的整数倍']); %h为步长（b-a）取为h的整数倍;
    y=0;
    for i=1:n-1
    	y=y+f(a+i*h);
    end
    y=y+0.5*(f(a)+f(b));%中间项求和
    fx=h*y;
end