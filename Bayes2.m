
%版权为刘文所有，2008年4月20日
%算法视线见《模式识别》P33－P44（各类样本的协方差不相等）
%为了提高实验样本测试的精度，故采用多次模拟求平均值的方法
N=input('实验模拟次数 N(N最好为奇数)＝ ');
Result(1:3,1:3)=0;      %判别矩阵的初始化
for k=1:N             %控制程序模拟次数N
    %生成二维正态分布的样本2 X N 维的矩阵 
    X1=mvnrnd([1 2],[4 0;0 6],300)';   %2 X N
    X2=mvnrnd([5 3],[5 0;0 1],200)';
    X3=mvnrnd([4 7],[2 0;0 9],500)';   %样本程序
    %---------------------------------------------------%
    %测试样本
    X10=mvnrnd([1 2],[4 0;0 6],100)';   %2 X N
    X20=mvnrnd([5 3],[5 0;0 1],100)';
    X30=mvnrnd([4 7],[2 0;0 9],100)';  
    %先验概率
    P(1)=length(X1)/(length(X1)+length(X2)+length(X3));
    P(2)=length(X2)/(length(X1)+length(X2)+length(X3));
    P(3)=length(X3)/(length(X1)+length(X2)+length(X3));
    %计算相关量  cov(X)：协方差矩阵 Ave:均值
    %--------------------------------------------------------%
    W1=-1/2*inv(cov(X1')); W2=-1/2*inv(cov(X2')); W3=-1/2*inv(cov(X3'));%
    Ave1=(sum(X1')/length(X1))';Ave2=(sum(X2')/length(X2))';
    Ave3=(sum(X3')/length(X3))';%计算平均值(2维列向量)
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
                    Result(1,2)=Result(1,2)+1;%记录误判情况
                else
                    Result(1,3)=Result(1,3)+1;%记录误判情况
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
%画出各样本的分布情况
subplot(2,1,1)
plot(X1(1,:),X1(2,:),'r.','LineWidth',2),hold on
plot(X2(1,:),X2(2,:),'go','LineWidth',2),hold on
plot(X3(1,:),X3(2,:),'b+','LineWidth',2),hold on
title('训练样本分布情况')
legend('训练样本1','训练样本2','训练样本3')
subplot(2,1,2)
plot(X10(1,:),X10(2,:),'r.','LineWidth',2),hold on
plot(X20(1,:),X20(2,:),'go','LineWidth',2),hold on
plot(X30(1,:),X30(2,:),'b+','LineWidth',2),hold on
title('测试样本分布情况')
legend('测试样本1','测试样本2','测试样本3')
%由于多次循环后存在小数，根据实际情况判别矩阵须取整
%如果N为偶数，可能出现小数为0.5的情况，此时将无法更加准确判断矩阵
Result=Result/N     %判别矩阵，反映Bayes的判别效果
for i=1:length(Result)
    if round(sum(Result(i,:)-fix(Result(i,:))))==1
        [m,n]=find(max(Result(i,:)-fix(Result(i,:)))==(Result(i,:)-fix(Result(i,:))));
        n=min(n);%存在小数点相同的情况随即选取一个
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
        n1=min(n1);n2=min(n2);%如果有存在小数点相同的情况,随即选取一个
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
