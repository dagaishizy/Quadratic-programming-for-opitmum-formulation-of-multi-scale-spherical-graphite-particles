function [fx] = Simpson(f,a,b,h)
%   复化辛普森公式
%   Input：f为要求导的函数；a,b为所求的函数的积分区间;h为步长,其中（b-a）取为h/2的整数倍;
%   Output：fx为f在[a,b]上的积分值
    assert(nargin == 4);
    n = round((b-a)/(2*h));
    %assert((rem(n,1) == 0),['更改步长h,使得(b-a)取为h/2的整数倍']); %h为步长（b-a）取为h/2的整数倍;
    fx = h/3 * (f(a) + 4*sum(f(a+h:2*h:a+(2*(n-1)*h))) + 2*sum(f(a+2*h:2*h:a+(2*(n-1)*h))) + f(b));
end