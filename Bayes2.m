
%��ȨΪ�������У�2008��4��20��
%�㷨���߼���ģʽʶ��P33��P44������������Э�����ȣ�
%Ϊ�����ʵ���������Եľ��ȣ��ʲ��ö��ģ����ƽ��ֵ�ķ���
N=input('ʵ��ģ����� N(N���Ϊ����)�� ');
Result(1:3,1:3)=0;      %�б����ĳ�ʼ��
for k=1:N             %���Ƴ���ģ�����N
    %���ɶ�ά��̬�ֲ�������2 X N ά�ľ��� 
    X1=mvnrnd([1 2],[4 0;0 6],300)';   %2 X N
    X2=mvnrnd([5 3],[5 0;0 1],200)';
    X3=mvnrnd([4 7],[2 0;0 9],500)';   %��������
    %---------------------------------------------------%
    %��������
    X10=mvnrnd([1 2],[4 0;0 6],100)';   %2 X N
    X20=mvnrnd([5 3],[5 0;0 1],100)';
    X30=mvnrnd([4 7],[2 0;0 9],100)';  
    %�������
    P(1)=length(X1)/(length(X1)+length(X2)+length(X3));
    P(2)=length(X2)/(length(X1)+length(X2)+length(X3));
    P(3)=length(X3)/(length(X1)+length(X2)+length(X3));
    %���������  cov(X)��Э������� Ave:��ֵ
    %--------------------------------------------------------%
    W1=-1/2*inv(cov(X1')); W2=-1/2*inv(cov(X2')); W3=-1/2*inv(cov(X3'));%
    Ave1=(sum(X1')/length(X1))';Ave2=(sum(X2')/length(X2))';
    Ave3=(sum(X3')/length(X3))';%����ƽ��ֵ(2ά������)
    w1=inv(cov(X1'))*Ave1;w2=inv(cov(X2'))*Ave2;w3=inv(cov(X3'))*Ave3;%2
    w10=-1/2*Ave1'*inv(cov(X1'))*Ave1-1/2*log(det(cov(X1')))+log(P(1));
    w20=-1/2*Ave2'*inv(cov(X2'))*Ave2-1/2*log(det(cov(X2')))+log(P(2));
    w30=-1/2*Ave3'*inv(cov(X3'))*Ave3-1/2*log(det(cov(X3')))+log(P(3));
    %-----------------------------------------------------------%
    for i=1:3                                     
        for j=1:100                               
            if i==1
                g1=X10(:,j)'*W1*X10(:,j)+w1'*X10(:,j)+w10;  
                g2=X10(:,j)'*W2*X10(:,j)+w2'*X10(:,j)+w20;
                g3=X10(:,j)'*W3*X10(:,j)+w3'*X10(:,j)+w30;
                if g1>=g2&g1>=g3     
                    Result(1,1)=Result(1,1)+1;
                elseif g2>=g1&g2>=g3
                    Result(1,2)=Result(1,2)+1;%��¼�������
                else
                    Result(1,3)=Result(1,3)+1;%��¼�������
                end
            elseif i==2
                g1=X20(:,j)'*W1*X20(:,j)+w1'*X20(:,j)+w10;
                g2=X20(:,j)'*W2*X20(:,j)+w2'*X20(:,j)+w20;
                g3=X20(:,j)'*W3*X20(:,j)+w3'*X20(:,j)+w30;
                if g2>=g1&g2>=g3
                    Result(2,2)=Result(2,2)+1;
                elseif g1>=g2&g1>=g3
                    Result(2,1)=Result(2,1)+1;
                else
                    Result(2,3)=Result(2,3)+1;
                end
            else
                g1=X30(:,j)'*W1*X30(:,j)+w1'*X30(:,j)+w10;
                g2=X30(:,j)'*W2*X30(:,j)+w2'*X30(:,j)+w20;
                g3=X30(:,j)'*W3*X30(:,j)+w3'*X30(:,j)+w30;
                if g3>=g1&g3>=g2
                    Result(3,3)=Result(3,3)+1;
                elseif g2>=g1&g2>=g3
                    Result(3,2)=Result(3,2)+1;
                else
                    Result(3,1)=Result(3,1)+1;
                end
            end
        end
    end
end
%�����������ķֲ����
subplot(2,1,1)
plot(X1(1,:),X1(2,:),'r.','LineWidth',2),hold on
plot(X2(1,:),X2(2,:),'go','LineWidth',2),hold on
plot(X3(1,:),X3(2,:),'b+','LineWidth',2),hold on
title('ѵ�������ֲ����')
legend('ѵ������1','ѵ������2','ѵ������3')
subplot(2,1,2)
plot(X10(1,:),X10(2,:),'r.','LineWidth',2),hold on
plot(X20(1,:),X20(2,:),'go','LineWidth',2),hold on
plot(X30(1,:),X30(2,:),'b+','LineWidth',2),hold on
title('���������ֲ����')
legend('��������1','��������2','��������3')
%���ڶ��ѭ�������С��������ʵ������б������ȡ��
%���NΪż�������ܳ���С��Ϊ0.5���������ʱ���޷�����׼ȷ�жϾ���
Result=Result/N     %�б���󣬷�ӳBayes���б�Ч��
for i=1:length(Result)
    if round(sum(Result(i,:)-fix(Result(i,:))))==1
        [m,n]=find(max(Result(i,:)-fix(Result(i,:)))==(Result(i,:)-fix(Result(i,:))));
        n=min(n);%����С������ͬ������漴ѡȡһ��
        for j=1:length(Result)
            if j==n
                Result(i,j)=fix(Result(i,j))+1;
            else
                Result(i,j)=fix(Result(i,j));
            end
        end
    elseif round(sum(Result(i,:)-fix(Result(i,:))))==2
        [m,n1]=find(max(Result(i,:)-fix(Result(i,:)))==(Result(i,:)-fix(Result(i,:))));
        [m,n2]=find(min(Result(i,:)-fix(Result(i,:)))==(Result(i,:)-fix(Result(i,:))));
        n1=min(n1);n2=min(n2);%����д���С������ͬ�����,�漴ѡȡһ��
        for j=1:length(Result)
            if j==n1
                Result(i,j)=fix(Result(i,j))+1;
            elseif j==n2
                Result(i,j)=fix(Result(i,j));
            else
                Result(i,j)=fix(Result(i,j))+1;
            end
        end
    else
        continue,
    end
end
