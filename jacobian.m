

syms Q1 Q3 H1 H2 H3 s
dh1 = 64.935*Q1 - (7.216*10^(-3))*(sqrt(H1 - H2));
dh2 = 7.216*10^(-3)*(sqrt(H1 - H2)+(sqrt(H3 - H2))) - 0.0101*sqrt(H2);
dh3 = 64.935*Q3 - (7.216*10^(-3))*(sqrt(H3 - H2));

jac_A = jacobian([ dh1, dh2, dh3],[ H1, H2, H3]);
B = jacobian([ dh1, dh2, dh3],[ Q1, Q3]);
A = subs(jac_A,{ H1, H2, H3},{ 0.3, 0.2014, 0.3});

B2 = B(:,1);
B1 = B(:,2);

vec_i = [1 1 1];
II = diag(vec_i);

C1 = [1 0 0];
C2 = [0 0 1];

HS1 = C1*((s*II - A )^(-1))*B1;
HS2 = C1*((s*II - A )^(-1))*B2;
HS3 = C2*((s*II - A )^(-1))*B1;
HS4 = C2*((s*II - A )^(-1))*B2;

HS11 = vpa(HS1, 3);
HS22 = vpa(HS2, 3);
HS33 = vpa(HS3, 3);
HS44 = vpa(HS4, 3);

pretty(HS11)
pretty(HS22)
pretty(HS33)
pretty(HS44)

