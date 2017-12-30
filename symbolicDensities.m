% Density Matrix Symbolic Calculations
%
% mikael.mieskolainen@cern.ch, 2017
clear; close all;


%% Spin-1 Analytical Case

a = sym('a');
b = sym('b');
c = sym('c');
d = sym('d');

% Construct spin-1 initial state density matrix with parity conservation
rho1 = [(1-a)/2, b+1i*c, d;
        b-1i*c,   a, -b+1i*c;
        d, -b-1i*c, (1-a)/2];

% Calculate traces s_k
s = cell(3,1);
for k = 1:length(s)
    s{k} = simplify(trace(rho1^k));
end

% Positivity conditions, >= 0
% P. Minnaert, Physical Review, 1966
pos_cond = cell(4,1);

pos_cond{1} = simplify( -s{2}+1 );
pos_cond{2} = simplify( 2*s{3}-3*s{2}+1 );

% Solve the eigenvalues
lambda1 = simplify(eig(rho1))


%% Spin-2 Analytical Case

a = sym('a');
b = sym('b');
c = sym('c');
d = sym('d');
e = sym('e');
f = sym('f');
g = sym('g');
h = sym('h');
j = sym('j');
k = sym('k');
l = sym('l');
m = sym('m');

% Construct spin-2 initial state density matrix with parity conservation
rho2 = [a,       f+1i*c,   g+1i*d,      h+1i*e,      m;
       f-1i*c,   b,        j+1i*k,      l,          -h+1i*e;
       g-1i*d,   j-1i*k,   1-2*(a+b),  -j+1i*k,      g-1i*d;   
       h-1i*e,   l,       -j-1i*k,      b,          -f+1i*c;
       m,       -h-1i*e,   g+1i*d,     -f-1i*c,      a];

% Calculate traces s_k
s = cell(5,1);
for k = 1:length(s)
    s{k} = simplify(trace(rho2^k));
end

% Positivity conditions,  >= 0
% P. Minnaert, Physical Review, 1966
pos_cond = cell(4,1);

pos_cond{1} = simplify( -s{2}+1 );
pos_cond{2} = simplify( 2*s{3}-3*s{2}+1 );
pos_cond{3} = simplify( -6*s{4}+8*s{3}+3*(s{2})^2-6*s{2}+1 );
pos_cond{4} = simplify( 24*s{5}-30*s{4}+20*s{3}-20*s{3}*s{2}+15*(s{2})^2-10*s{2}+1 );


% Try to solve the eigenvalues
% -> too high dimensional polynomial (over 4) for a closed form solution
lambda2 = simplify(eig(rho2));


%% Spin-2 density matrix

% Basis vectors
B = eye(5);

% Quantum superposition state
cc = [0.25+0.1i 0.1+0.05i 0.3+0.1i 0.1+0.05i 0.25+0.002i];

q = zeros(5,1);
for i = 1:5
    q = q + B(:,i)*cc(i);
end
% Normalization
q = q / norm(q,2);

% Construct density matrix
rho = q*q'

