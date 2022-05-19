function [fx] = Simpson(f,a,b,h)
%   ��������ɭ��ʽ
%   Input��fΪҪ�󵼵ĺ�����a,bΪ����ĺ����Ļ�������;hΪ����,���У�b-a��ȡΪh/2��������;
%   Output��fxΪf��[a,b]�ϵĻ���ֵ
    assert(nargin == 4);
    n = round((b-a)/(2*h));
    %assert((rem(n,1) == 0),['���Ĳ���h,ʹ��(b-a)ȡΪh/2��������']); %hΪ������b-a��ȡΪh/2��������;
    fx = h/3 * (f(a) + 4*sum(f(a+h:2*h:a+(2*(n-1)*h))) + 2*sum(f(a+2*h:2*h:a+(2*(n-1)*h))) + f(b));
end