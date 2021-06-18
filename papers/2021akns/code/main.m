clear;clc;
%% Inital potential functions
%Example 1
p = @(x) sin(4*pi.*x);
q = @(x) cos(2*pi.*x);

%% Inital hyperparameter
N = 30;                 % Determine the number of half of spectra E
K = 60;                 % Determine the number of zeros of Y_plus,Z_plus
P = fix(K/3)-1;         % Determine the number of zeros of Y0,Z0
M = P-1;                % Determine the number of Drichlet eigenvalues and Norming constants

%% Find half of spectra data E of (p,q) on [0,1], length(E)=2N+1.
model_E = AKNS_odemodel([0,1],p,q);
E_xspan = [-N,N]*2*pi;
while model_E.solution_zeros(E_xspan(1)) < 2*N+2
    E_xspan(1) = E_xspan(1) - 1;
end
while model_E.solution_zeros(E_xspan(2)) < 2*N+2
    E_xspan(2) = E_xspan(2) + 1;
end
while model_E.solution_zeros(E_xspan(1)) >= 2*N+2
    E_xspan(1) = E_xspan(1) + 1;
end
while model_E.solution_zeros(E_xspan(2)) >= 2*N+2
    E_xspan(2) = E_xspan(2) - 1;
end
spectra = find_zeros(model_E.Y,E_xspan);
E = spectra(1:2:end);
clear spectra E_xspan model_E

%% Solve F+(1/4,e), then find zeros A,B of Y+(1/4,e), Z+(1/4,e),length(A)=2K+1,length(B)=2K
model_AB = AKNS_odemodel([1,1/4],p,q);
A_xspan = [-K,K]*(4/3)*pi;
while model_AB.solution_zeros(A_xspan(1)) < K+2
    A_xspan(1) = A_xspan(1) - 1;
end
while model_AB.solution_zeros(A_xspan(2)) < K+2
    A_xspan(2) = A_xspan(2) + 1;
end
while model_AB.solution_zeros(A_xspan(1)) >= K+2
    A_xspan(1) = A_xspan(1) + 1;
end
while model_AB.solution_zeros(A_xspan(2)) >= K+2
    A_xspan(2) = A_xspan(2) - 1;
end
A = find_zeros(model_AB.Y,A_xspan);

B_xspan = [-K,K-1]*(4/3)*pi+2/3*pi;
while model_AB.solution_zeros(B_xspan(1)) < K+1
    B_xspan(1) = B_xspan(1) - 1;
end
while model_AB.solution_zeros(B_xspan(2)) < K+1
    B_xspan(2) = B_xspan(2) + 1;
end
while model_AB.solution_zeros(B_xspan(1)) >= K+1
    B_xspan(1) = B_xspan(1) + 1;
end
while model_AB.solution_zeros(B_xspan(2)) >= K+1
    B_xspan(2) = B_xspan(2) - 1;
end
B = find_zeros(model_AB.Z,B_xspan);
Y_plus = @(e) model_AB.Y(e);
Z_plus = @(e) model_AB.Z(e);
clear A_xspan B_xspan model_AB

%% Solve W1(e)
coeff_W1 = sin(N*2*pi)/gen_W1(E,N*4*pi);
W1 = @(e) coeff_W1*gen_W1(E,e);
clear coeff_W1

%% Solve Y0,Z0 and their zeros C,D
diff = @(f,e) (f(e+0.01)-f(e-0.01))/0.02;

remainder_Y0 = -W1(A)./Z_plus(A) - sin(1/4*A);
diff_Y_plus_on_A = diff(Y_plus,A);
fun_Y0 = @(e) sin(1/4*e) + sum(Y_plus(e)./(e-A)./diff_Y_plus_on_A.*remainder_Y0);
interpolation_point_Y0 = linspace(A(1)-30,A(end)+30,1000);
value_interpolation_point_Y0 = fun_Y0(interpolation_point_Y0);
Y0 = @(e) (e<interpolation_point_Y0(1) | e>interpolation_point_Y0(end)).*sin(1/4*e) + (e>=interpolation_point_Y0(1) & e<=interpolation_point_Y0(end)).*interp1(interpolation_point_Y0,value_interpolation_point_Y0,e,'linear',0);

remainder_Z0 = W1(B)./Y_plus(B) - cos(1/4*B);
diff_Z_plus_on_B = diff(Z_plus,B);
fun_Z0 = @(e) cos(1/4*e) + sum(Z_plus(e)./(e-B)./diff_Z_plus_on_B.*remainder_Z0);
interpolation_point_Z0 = linspace(B(1)-30,B(end)+30,1000);
value_interpolation_point_Z0 = fun_Z0(interpolation_point_Z0);
Z0 = @(e) (e<interpolation_point_Z0(1) | e>interpolation_point_Z0(end)).*cos(1/4*e) + (e>=interpolation_point_Z0(1) & e<=interpolation_point_Z0(end)).*interp1(interpolation_point_Z0,value_interpolation_point_Z0,e,'linear',0);

CD_xspan = [-P*4*pi,P*4*pi];
C = find_zeros(Y0,CD_xspan);
D = find_zeros(Z0,CD_xspan);
clear remainder_Y0 diff_Y_plus_on_A fun_Y0 remainder_Z0 diff_Z_plus_on_B fun_Z0 CD_xspan ...
interpolation_point_Y0 value_interpolation_point_Y0 interpolation_point_Z0 value_interpolation_point_Z0

%% Solve F-(1/4,e)
remainder_Y_minus = -Y_plus(C) - sin(1/4*C);
diff_Y0_on_C = diff(Y0,C);
fun_Y_minus = @(e) sin(1/4*e) + sum(Y0(e)./(e-C)./diff_Y0_on_C.*remainder_Y_minus);
interpolation_point_Y_minus = linspace(C(1)-10,C(end)+10,1000);
value_interpolation_point_Y_minus = fun_Y_minus(interpolation_point_Y_minus);
Y_minus = @(e) (e<interpolation_point_Y_minus(1) | e>interpolation_point_Y_minus(end)).*sin(1/4*e) + (e>=interpolation_point_Y_minus(1) & e<=interpolation_point_Y_minus(end)).*interp1(interpolation_point_Y_minus,value_interpolation_point_Y_minus,e,'linear',0);

remainder_Z_minus = -Z_plus(D) - cos(1/4*D);
diff_Z0_on_D = diff(Z0,D);
fun_Z_minus = @(e) cos(1/4*e) + sum(Z0(e)./(e-D)./diff_Z0_on_D.*remainder_Z_minus);
interpolation_point_Z_minus = linspace(D(1)-10,D(end)+10,1000);
value_interpolation_point_Z_minus = fun_Z_minus(interpolation_point_Z_minus);
Z_minus = @(e) (e<interpolation_point_Z_minus(1) | e>interpolation_point_Z_minus(end)).*cos(1/4*e) + (e>=interpolation_point_Z_minus(1) & e<=interpolation_point_Z_minus(end)).*interp1(interpolation_point_Z_minus,value_interpolation_point_Z_minus,e,'linear',0);

%% check F-(1/4,e)
model = AKNS_odemodel([0,1/4],p,q);
exact_Y = @(e) model.Y(e);
exact_Z = @(e) model.Z(e);
x = linspace(C(1)/2,C(end)/2,200);
y1 = exact_Y(x);
y2 = Y_minus(x);
z1 = exact_Z(x);
z2 = Z_minus(x);
figure(1)
plot(x,y1,'r',x(1:2:end),y2(1:2:end),'ob')
title('The function Y_{-}(1/4,\lambda) plot')
xlabel('\lambda')
ylabel('Y_{-}(1/4,\lambda)')
legend('Exact','Numerical')

figure(2)
plot(x,z1,'r',x(1:2:end),z2(1:2:end),'ob')
title('The function Z_{-}(1/4,\lambda) plot')
xlabel('\lambda')
ylabel('Z_{-}(1/4,\lambda)')
legend('Exact','Numerical')

%% Fit (p,q) on the [0,1/4]
fun_E = @(e) exp(1/4*1i*e).*(Z_minus(e)-1i*Y_minus(e));
fun_S = @(e) conj(fun_E(e))./fun_E(e);

e1 = -4*pi-2*rand(1);
while abs(Y_minus(e1)-sin(1/4*e1)) + abs(Z_minus(e1)-cos(1/4*e1)) > 1e-5
    e1 = e1 - 4*pi;
end
e2 = 4*pi+5*rand(1);
while abs(Y_minus(e2)-sin(1/4*e2)) + abs(Z_minus(e2)-cos(1/4*e2)) > 1e-5
    e2 = e2 + 4*pi;
end

fun_F = @(x) 1/(2*pi)*integral(@(e) (1-fun_S(e)).*exp(1i*e.*x),e1,e2);

N = 200;
x = linspace(0,1/4,N+1);
h = x(2)-x(1);

F_value = zeros(2*N+1,1);
for i = 1:2*N+1
    F_value(i) = fun_F((i-1)*h);
end
F = zeros(N+1,N+1);
for i = 1:N+1
    F(i,:) = F_value(i:N+i);
end

F_A1 = triu(F);
F_A2 = F_A1 - 1/2*diag(diag(F_A1));
F_A2(:,N+1) = F_A2(:,N+1)./2;
F_A2(N+1,N+1) = 0;
F_A = h*F_A2;

K_f1 = zeros(N+1,N+1);
K_f1_new = -conj(F)-conj(F_A*K_f1);
while norm(K_f1_new-K_f1) > 1e-6
    K_f1 = K_f1_new;
    K_f1_new = -conj(F)-conj(F_A*K_f1);
end
K11 = real(diag(K_f1));
K12 = imag(diag(K_f1));

K_f2 = zeros(N+1,N+1);
K_f2_new = -conj(1i*F)-conj(F_A*K_f2);
while norm(K_f2_new-K_f2) > 1e-6
    K_f2 = K_f2_new;
    K_f2_new = -conj(1i*F)-conj(F_A*K_f2);
end
K21 = real(diag(K_f2));
K22 = imag(diag(K_f2));

y_p = K11-K22;
y_q = K12+K21;

figure(3)
plot(x,p(x),'linewidth',1.1,'color',[252,41,30]/255)
hold on
plot(x,y_p,'o','MarkerIndices',1:5:length(x),'Markersize',5,'markerfacecolor', [0, 70, 222]/255)
legend('Exact p','Numerical p')

figure(4)
plot(x,q(x),'linewidth',1.1,'color',[252,41,30]/255)
hold on 
plot(x,y_q,'o','MarkerIndices',1:5:length(x),'Markersize',5,'markerfacecolor', [0, 70, 222]/255)
legend('Exact q','Numerical q')

%% plot figure2 and save plot data
%print(1,'-depsc','-r300','function_y');
%print(2,'-depsc','-r300','function_z');
print(3,'-depsc','-r300','p');
print(4,'-depsc','-r300','q');
save plotData.mat x p q y_p y_q
