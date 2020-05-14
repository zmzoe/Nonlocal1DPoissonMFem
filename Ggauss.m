function G=Ggauss(h,g)
GaussP=[-0.7745967 0 0.7745967];                            %��˹��
GaussA=[0.5555556 0.8888889 0.5555556];                     %��˹ϵ��

XX = 0:h:1;                                                  %����[0,1]
n = length(XX)-1;
G=zeros(2*n-1,1);

%% odd  �������� g*(phi_2i-1 - phi_2i)
for i=1:n
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %����任(Ѱ��[x_i,x_i1]����Ӧ�Ļ��ֵ�)
    
    
    w =@(x) g(x)*(XX(i+1)+XX(i)-2*x)/h;                              %�������� g*(phi_2i-1 - phi_2i)
     
    for k=1:3
        G(2*i-1) = G(2*i-1) + h/2*w(points(k))*GaussA(k);
    end
end

%% even �������� g*(phi_2i)��[x_i,x_i1], ��ȥ g*(phi_2i1)��[x_i1,x_i2]

% part 1-------�������� g*(phi_2i)��[x_i,x_i1]

G1 = zeros(n-1,1);
G2 = zeros(n-1,1);
for i=1:n-1
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %����任(Ѱ��[x_i,x_i1]����Ӧ�Ļ��ֵ�)
    
   
    w =@(x) g(x)*(x-XX(i))/h;                                       % �������� g*(phi_2i)
     
    for k=1:3
        G1(i) = G1(i) + h/2*w(points(k))*GaussA(k);
    end
end

% part 2-------�������� g*(phi_2i1)��[x_i1,x_i2]
for i=1:n-1
    points = h/2*GaussP + (XX(i+2)+XX(i+1))/2;                  %����任(Ѱ��[x_i,x_i1]����Ӧ�Ļ��ֵ�)
    
   
    w =@(x) g(x)*(XX(i+2)-x)/h;                                       %�������� g*(phi_2i1)
    
    for k=1:3
        G2(i) = G2(i) + h/2*w(points(k))*GaussA(k);
    end
end

for i=1:n-1
    G(2*i)= G1(i)-G2(i);
end

