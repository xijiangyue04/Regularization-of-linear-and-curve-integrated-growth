%���������input_pnts(nx3)  ����ʽ������TM  ����κ���X=a(1)T^2+b(1)T+c(1);Y=a(2)T^2+b(2)T+c(2),... TM=2
%���������parameter(3xTM) ��һ��ΪX������ϲ������ڶ���ΪY������ϲ�����������ΪZ������ϲ���; ��������ϵĵ�curvefit_pnts(nx3)
%mean_dis ��Ӧ�����ƫ��ƽ��ֵ
function [parameter,mean_dis] = space_curve_LS(input_pnts,TM)

for i=1:4
    if i<4
        pnts=sortrows(input_pnts,i,'descend');
    end
    if i==4
        pnts=input_pnts;
    end
X=pnts(:,1)';
Y=pnts(:,2)';
Z=pnts(:,3)';
%������ɢ������ά����X,Y,Z��ͬ��
n=size(pnts,1);
T=0:1:(n-1)*1;
%��X,Y,Z�ֱ���T֮����ж���ʽ��ϣ��������ʽ�����Ϊ8��
C1=polyfit(T,X,TM);  %X����������ϲ��� ϵ������������,��X=C1(1)*T^2+C1(2)*T+C1(3)
C2=polyfit(T,Y,TM);  %Y����������ϲ���
C3=polyfit(T,Z,TM);   %Z����������ϲ���
ZX=polyval(C1,T); %������Ϻ�Ĺ���x��������(1xn)
ZY=polyval(C2,T);%������Ϻ�Ĺ���y��������(1xn)
ZZ=polyval(C3,T);%������Ϻ�Ĺ���z��������(1xn)
curvefit_pnts=[ZX;ZY;ZZ]';
fit_dis=sqrt(sum((curvefit_pnts-pnts).^2,2)); %nx1
mean_dis(i,:)=mean(fit_dis);
parameter_total{i}=[C1;C2;C3];
end
parameter=parameter_total{mean_dis==min(mean_dis)};






