clear all;
clc;
format long;
syms u v 

r = 4.0;
No = 30;
d = 1.3;
base = 1.0;

f =  u * ((1.0 - u - v) * ((r * u / (u + v)) * (1.0 - ((1.0 - (1.0 - (1.0 - u - v))^No) / (No * (u + v)))) ...
    - (1.0 + (r - 1.0) * (1.0 - u - v)^(No - 1) - (r / No) * ((1.0 - (1.0 - u - v)^No) / (1.0 - (1.0 - u - v)))) + base) - d); %% activator_eq
g = v * ((1.0 - u - v) * ((r * u / (u + v)) * (1.0 - ((1.0 - (1.0 - u - v)^No) / (No * (u + v)))) + base) - d); %%% inactivator_eq

ue = 0.020078; 
ve = 0.102254; %% equilibrium





% Derivatives of f
fu = diff(f, u);
fv = diff(f, v);

% Derivatives of g
gu = diff(g, u);
gv = diff(g, v);  

%%%%% substitution of equilibrium 
fu_sub = subs(fu, [u, v], [ue, ve]); 
fv_sub = subs(fv, [u, v], [ue, ve]);
gu_sub = subs(gu, [u, v], [ue, ve]);
gv_sub = subs(gv, [u, v],[ue, ve]);

%%%%%%% converting syms to double 
fu_num = double(fu_sub);
fv_num = double(fv_sub);
gu_num = double(gu_sub);
gv_num = double(gv_sub);


epsilon = 0.01; %%%%%% Diffusion of activator
L=load("laplacian_scalefree_N=1000,k=10.dat"); %% load_laplacian
e=eig(L); %% laplacian eigen value
sigma = 25;
data2=[];

for j=1:length(e)
     Lam =-e(j); %%% making eigenvalues negative since python and other software gives eigenvalues of the laplacian positive
 
    %%%%%%% calculation of dispersion
    % First term
trace1 = (fu_num + gv_num + (1 + sigma) * epsilon * Lam);
deter1 = (fu_num * gv_num - fv_num * gu_num) + epsilon * Lam * ( sigma * fu_num + gv_num) + sigma * epsilon ^2 * Lam^2;

% Square root term
%term2 = sqrt(4 * fv_num * gu_num + (fu_num - gv_num + (1 - sigma) * epsilon * Lam)^2);


lambda_plus = (trace1 + sqrt(trace1^2 - 4*deter1))/2;
%%%%%%%%%%%%%%%
%lambda_plus = (fu_num * gv_num - fv_num * gu_num) + epsilon * Lam * ( sigma * fu_num + gv_num) + sigma * epsilon ^2 * Lam^2;

data2=[data2; Lam real(lambda_plus)]; %%%%%%% storing the dispersion rate for each eigenvalues of the Laplacian
    end

figure;
plot((-data2(:,1)), data2(:,2))
hold on;
%plot((-data2(:,1)), data2(:,2),'ro')
%yline(0)

%hold on;
%plot([0,30],[0,0])


