function [fx]=compositeTrapezoidintegrationrule(f,a,b,h)
%   �������ι�ʽ
%   Input��fΪҪ�󵼵ĺ�����a,bΪ����ĺ����Ļ�������;hΪ����,���У�b-a��ȡΪh��������;
%   Output��fxΪf��[a,b]�ϵĻ���ֵ
    assert(nargin==4);
    n = round((b-a)/h);
%assert((rem(n,1)==0),['���Ĳ���h,ʹ��(b-a)ȡΪh��������']); %hΪ������b-a��ȡΪh��������;
    y=0;
    for i=1:n-1
    	y=y+f(a+i*h);
    end
    y=y+0.5*(f(a)+f(b));%�м������
    fx=h*y;
end