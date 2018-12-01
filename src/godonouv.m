%% 線形双曲線型連立偏微分方程式　ゴドノフ法 2x2
% dU/dt + dF/dx = 0, F = AU
% U = 2x1, A = 2x2, F = 2x1
close all;

A = [3,1;2,4]; 
n= 100;
x = zeros([1,n]);
U = zeros([2,n]);
U_init = zeros([2,n]);
Flux = zeros([2,n]); %Flux(j) = F(j+1/2)

nlast = 100; %number of timesteps
cfl = 0.5;
dx = 1/n;
lambda = eig(A);
a= max(lambda);
dt = cfl * dx / a;
disp('dt:');
disp(dt);
%% Initialize
U_init_l = [1;1];
U_init_r = [2;-1];

for i = 1:n
    x(i) = (i-1)*dx;
    if(x(i) <= 0.2)
        U(:,i) = U_init_l;
        U_init(:,i) = U_init_l;
    else
        U(:,i) = U_init_r;
        U_init(:,i) = U_init_r;
    end
end
%% Exact Answer
[R,Lambda] = eig(A);
disp('Lambda:');
disp(Lambda);
disp('R:');
disp(R);
W_init_l = R\U_init_l;
W_init_r = R\U_init_r;
W_exact = zeros([2,100]);

for j = 1:n
    for k = 1:2
       if(x(j) < (0.2 + dt*lambda(k)*nlast))
           W_exact(k,j) = W_init_l(k,1);
       else
           W_exact(k,j) = W_init_r(k,1);
       end
    end
end
U_exact = R * W_exact;
%% Numerical Answer
for t = 1:nlast
    F = A * U;
    U_new = zeros([2,n]); %U of next timestep
    for m = 1:n-1
        Flux(:,m) = 0.5 * (F(:,m) + F(:,m+1) - R * (abs(Lambda) * inv(R)) * (U(:,m+1) - U(:,m)));
    end
    for p = 2:n-1
        U_new(:,p) = U(:,p) - (dt/dx)*(Flux(:,p) - Flux(:,p-1));
    end
    U_new(:,1) = repmat(U(:,1),1);
    %U_new(:,2) = repmat(U(:,2),1);
    U_new(:,n) = repmat(U(:,n),1);
    U = repmat(U_new,1);
end      
%% Plot Answer
figure(1)
plot(x,U_init(1,:),'r--'); hold on;
plot(x,U_exact(1,:),'r'); hold on;
plot(x,U(1,:),'r-o'); hold on;
plot(x,U_init(2,:),'b--'); hold on;
plot(x,U_exact(2,:),'b'); hold on;
plot(x,U(2,:),'b-o'); hold on;

legend('u(initial)','u(exact)','u(numerical)','v(initial)','v(exact)','v(numerical)','Location','best');
title('CFL=0.5,100steps');
xlabel('x');
ylabel('u,v');
xlim([0,1]);
ylim([-1.5,2.5]);