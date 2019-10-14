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

sys = ss(A,B,C,D);
e = eig(sys);

dt = 0.1;
t = [0:dt:500];
u = [zeros(size(t)) ; ones(size(t))];

zeta = abs(real(e))./abs(e);
phi = acos(zeta);
wn = abs(e).*cos(phi)./zeta;
wa = wn.*sqrt(1-zeta.^2);
tp = pi./wa;
Mp = 100.*exp(-pi./tan(phi));
ts = 4./(zeta.*wn);
disp('Stablilité des modes')
disp('    Zeta      Wn        Wa        tp        Mp        ts      ')
disp([zeta wn wa tp Mp ts])

%% paramètres du root locus
ft = tf(sys);

ft_a_v = ft(1,2); %% v/a_prop
z = zero(ft_a_v);
p = pole(ft_a_v);
sigma = sum(p)-sum(z);
[num den] = tfdata(ft_a_v, 'v');
syms s;
nums = poly2sym(num,s);
numd = diff(nums, s);
dens = poly2sym(den,s);
dend = diff(dens, s);
eq1 = -(nums*dend-dens*numd)/(nums^2);
breakin = vpa(solve(eq1,s), 3);
for n = 1:numel(p) 
    depart(n) = wrapToPi(pi - sum(angle(p(n)-p)) + sum(angle(p(n)-z)));
    arrivee(n) = wrapToPi(pi + sum(angle(p(n)-p)) - sum(angle(p(n)-z)));
end

angles = [num2str(rad2deg(depart)') repmat('  ',[4 1]) num2str(rad2deg(arrivee)')];
polestr = num2str(p);

disp('Paramètres du root locus')
disp('   Pole              départ         arrivée')
disp([polestr repmat('  ',[4 1]) angles])
disp([newline '   Break-in = ' char(vpa(breakin(1),3))])

k = [0:0.001:100];

figure
rlocus(ft(1,2),k)
%% 
Kv = 0.01
B(:,2)*Kv*C(1,:)
A1 = A + B(:,2)*Kv*C(1,:)
B1 = B(:,1)
C1 = C
D1 = D(:,1)
sys1 = ss(A1,B1,C1,D1)
ft_a_v2 = feedback(ft_a_v, Kv)

figure
step(sys, t)
figure
step(sys1, t)

%% marges de gain et phase

margin(ft_a_v2)


%% Réduction d'ordre

[num den] = tfdata(ft_a_v, 'v');
[r p k] = residue(num, den);

abs(r)./real(p)

[num1 den1] = residue(r(3:4), p(3:4), k)
[num2 den2] = residue(r(1:2), p(1:2), k)
H1 = tf(num1, den1)
H2 = tf(num2, den2)
figure
step(ft_a_v, H1, H2, t)

Kv = 0.1
num_a = [num1(1) num1(2)]
den_a = [den1(1) den1(2)+num1(1).*Kv den1(3)+num1(2).*Kv]

H_a = tf(num_a, den_a)
figure(7)
step(H_a, t)

figure
rlocus(H1)

%% Boucle externe

H = tf(sys1);
H5 = H(5)
rlocus(H(5))

[num den] = tfdata(H5, 'v')
syms w s Kp
symH5 = poly2sym(num, s)/poly2sym(den, s)
symH5 = subs(symH5, s, i*w)
eq = 1 + Kp*symH5 == 0
real(eq)



