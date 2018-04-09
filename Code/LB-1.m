f = fopen('C:\input.txt'); %input.txt

for i = 1:5 %считываем информацию с файла
    
    str = '';
    C = textscan(f, '%s', 1);
   
    for k = 1:length(C)
        str = [str char(C{k}) ' '];
    end;
        eval(['mass',num2str(i), ' = str2num(str);']);
end;

X = mass1';
U = mass2;
N = mass3;
Y = mass4;
T = mass5;

x = length(X);
u = length(U);
n = length(N);
y = length(Y);

A = randn(x, x);
B = randn(x, u);
E = randn(x, n);

C = randn(y, x);
D = randn(y, u);
F = randn(y, n);

t=0;
for t=1:T
    
    N = randn(n, 1);
    X = A*X + B*U' + E*N;
    Y = C*X + D*U' + F*N;

    Xtransf(:,t) = X;
    Ytransf(:,t) = Y;

end

dlmwrite('output.txt', Y);

%____проверка на устойчивость

onematrix = eye(x,x);
symbol = sym('k');
opredelitel = det(symbol*onematrix-A);
roots = solve(opredelitel,symbol);
roots = double(roots);

R=real(roots);
Im=imag(roots); 

flag1=0;
flag2=0;
flag3=0;

for i=1:length(roots)
if abs(R(i))<1
    flag1 = flag1+1;
else
    flag3=flag3+1;
end;
end

for i=1:length(Im)
if abs(Im(i))<1
    flag2=flag2+1;
else 
    flag3=flag3+1;
end;
end

if (flag1>0 || flag2>0) && (flag3<1)
    disp('Система устойчива')
else
    disp('Система не устойчива')
end;

%___время переходного процесса

minrealroots = min(R);
delta = 0.05;
tp = (1/abs(minrealroots))*(log(1/delta))

%___перерегулирование

Ymax = max(Y);
Yust = Y(length(Y),1);
Sigma = ((Ymax-Yust)/Yust)*100

%___фильтр

window = 3;
alpha=2/(window+1);

Ny = length(Ytransf);
Ny1 = length(Ytransf(:,1));
    for i=1:Ny
        Y1=zeros(Ny1,1);
        Y1=Ytransf(:,i);
        
        Filtery=zeros(1,Ny1);
        Filtery(1)=Y1(1);
        
        for j=2:Ny1
            Filtery(j) = alpha*Y1(j) + (1-alpha)*Filtery(j-1);
            Filter=Filtery'
        end
        Filterytr(:,i)=Filter;
    end

plot(1:T,Ytransf);
hold on;
plot(1:T, Filterytr);
