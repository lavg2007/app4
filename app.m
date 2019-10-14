clc
clear
close all
A = [-0.018223 -0.088571 -9.78 0 ;
     -0.003038 -1.2563 0 1;
     0 0 0 1;
     0.0617 -28.075 0 -4.5937];
B = [0 1.1962;
      0 -0.00120;
      0 0;
      7.84 -4.05];
C = [1 0 0 0;
       0 57.296 0 0;
       0 0 57.296 0;
       0 0 0 57.296;
       0 -57.296 57.296 0];
   
D = zeros(5,2);

dt = 0.1;
t = [0:dt:500];
u = [zeros(size(t)) ; ones(size(t))];

sys = ss(A,B,C,D);

[wn zeta p] = damp(sys);
phi = acos(zeta);
wa = wn.*sqrt(1-zeta.^2);
tp = pi./wa;
Mp = 100.*exp(-pi./tan(phi));
ts = 4./(zeta.*wn);


disp('Stablilité des modes')
disp('   Zeta      Wn        Wa        tp        Mp        ts      ')
disp([zeta wn wa tp Mp ts])

%% paramètres du root locus
ft = tf(sys);

ft_a_v = ft(1,2); %% v/a_prop

z = zero(ft_a_v);
p = pole(ft_a_v);
n = numel(p);
m = numel(z);

gamma = unique(wrapToPi((1+2.*[0:10]).*pi./(n-m)));
sigma = sum(p)-sum(z)/(n-m);
[num den] = tfdata(ft_a_v, 'v');
syms s;

nums = poly2sym(num,s);
numd = diff(nums, s);
dens = poly2sym(den,s);
dend = diff(dens, s);
eq1 = -(nums*dend-dens*numd)/(nums^2);

breakin = solve(eq1,s);
for n = 1:numel(p) 
    depart(n) = rad2deg(wrapToPi(pi - sum(angle(p(n)-p)) + sum(angle(p(n)-z))));
end

for n = 1:numel(z) 
    arrivee(n) = rad2deg(wrapToPi(pi + sum(angle(z(n)-p)) - sum(angle(z(n)-z))));
end

angles_p = num2str(depart');
polestr = num2str(p);
angles_z = num2str(arrivee');
zerostr = num2str(z);

disp(['Paramètres du root locus' newline])
disp(['   Gamma = ' num2str(gamma) newline]) 
disp(['   Sigma = ' num2str(sigma) newline]) 
disp(['   Break-in = ' char(vpa(breakin(2),3)) newline])
disp('   Pole              départ     ')
disp([repmat('  ',[4 1]) polestr repmat('  ',[4 1]) angles_p])
disp(newline)
disp('   Zéros             Arrivée')
disp([repmat('  ',[3 1]) zerostr repmat('  ',[3 1]) angles_z])


figure
rlocus(ft(1,2))
%% 
Kv = 1.03;
A1 = A - B(:,2)*Kv*C(1,:);
B1 = B(:,1);
C1 = C;
D1 = D(:,1);
sys1 = ss(A1,B1,C1,D1);
disp('Matrice A après ajout de boucle interne')
disp(A)


disp('Fonction de transfert V/a_prop, avec boucle interne Kv = 1.03')
ft_a_v1 = feedback(ft_a_v, Kv)

t = [0:0.01:10];

figure
step(ft_a_v1, t)


%% marges de gain et phase

%margin(ft_a_v2)


%% Réduction d'ordre

[num den] = tfdata(ft_a_v, 'v');
[r p k] = residue(num, den);

t = [0:0.01:100];

[num1 den1] = residue(r(3:4), p(3:4), k)
[num2 den2] = residue(r(1:2), p(1:2), k)
H1 = tf(num1, den1)
H2 = tf(num2, den2)
%figure
%step(ft_a_v, H1, H2, t)

[num den] = tfdata(H1, 'v');
syms s;

nums = poly2sym(num,s);
numd = diff(nums, s);
dens = poly2sym(den,s);
dend = diff(dens, s);
eq1 = -(nums*dend-dens*numd)/(nums^2);
breakin = solve(eq1,s);

syms sKv1
eq2 = 1 + sKv1*nums/dens == 0
eq2 = subs(eq2, s, breakin(1))
Kv1 = eval(solve(eq2, sKv1))

H_a = feedback(tf(num1, den1),Kv1)
%figure(7)
%step(H_a, t)

figure
hold on
rlocus(H_a)
%rlocus(Kv*ft_a_v)
%xlim([-1.75 -0.5])
%ylim([-1 1])

%% Boucle externe

H = tf(sys1);
H5 = H(5)
rlocus(H(5))

[num den] = tfdata(H5, 'v')
syms w s Kp
symH5 = poly2sym(num, s)/poly2sym(den, s)
symH5 = subs(symH5, s, i*w)
eq = vpa(1 + Kp*symH5 == 0,4)
eq1 = w^4 - 41.16*w^2 +52.51 + Kp*697.6 == 0
eq2 = -7.1*w^3 + 84.09*w + Kp*564.3*w == 0
rep = solve([eq1 eq2])
vpa(rep.w, 10)
eq3 = subs(eq2, w, 5.82)
K_limite = eval(solve(eq3))

K_marge = 10^(-18/20)

figure
rlocus(H(5))
figure
margin(H(5))
figure
margin(K_marge*H(5))






