f = fopen('C:\input.txt');
for i = 1:5 %считываем информацию с файла
    str = '';
    C = textscan(f, '%s', 1);
    for k = 1:length(C)
        str = [str char(C{k}) ' '];
    end;
    eval(['mass',num2str(i), ' = str2num(str);']);
end;
X = mass1;
U = mass2;
N = mass3;
Y = mass4;
T = mass5;

x = length(X);
u = length(U);
n = length(N);

A = randn(x);
B = randn(u, x);
E = randn(n, x);

t = 0;
while t < T
    X = X*A + U*B + N*E;
    t = t + 1;
end

C = randn(x);
D = randn(u, x);
F = randn(n, x);

Y = X*C + U*D + N*F;
dlmwrite('output.txt', Y);


onematrix = eye(6,6);
symbol = sym('k');
opredelitel = det(symbol*onematrix-A);
roots = solve(opredelitel,symbol);
roots = double(roots);
roots1 = abs(roots);

flag1=0;
flag2=0;

R=real(roots); %здесь менять имя
Im=imag(roots); %здесь менять имя

for i=1:length(R)
if ((R(i)>-1) && (R(i)<1))
% flag1 = flag1+1;
else
flag1=flag1+1;
end;
end
%-----------------------------------------
if Im == 0
    flag2=1;
end;

for i=1:length(Im)
if ((Im(i)>-1) && (Im(i)<1))
else
flag2=flag2+1;
end;
end

if flag1>0 && flag2>0
    disp('Система не устойчива')
else
    disp('Система устойчива')
end;

minrealroots = min(R);
delta = 0.05;
tp = (1/abs(minrealroots))*(log(1/delta))

Ymax = max(Y);
Yust = Y(1,length(Y));
Sigma = ((Ymax-Yust)/Yust)*100


window = 3;
alpha=2/(window+1);

Nx = length(X);
EMAx = zeros(1, Nx);
EMAx(1) = X(1);
for i=2:Nx
     EMAx(i) = alpha*X(i) + (1-alpha)*EMAx(i-1);
end

YwEMA = EMAx*C + U*D + N*F;
Ny = length(YwEMA);
EMAy = zeros(1, Ny);
EMAy(1) = YwEMA(1);
for i=2:Ny
     EMAy(i) = alpha*YwEMA(i) + (1-alpha)*EMAy(i-1);
end

plot(YwEMA);
hold on;
plot(EMAy);
