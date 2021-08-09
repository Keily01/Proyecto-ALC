# Proyecto de ALC
# El siguiente codigo resuelve el sistema lineal Ax=b
# para ello toma dos casos uno con matrices aleatorias y otro con matrices
# definidas positiva y simetrica, en ambos casos se generan matrices de 
# dimension 2 a 10, 20 a 100 (de 10 en 10) y de 200 a 900 (de 100 en 100).
# Para la resolucion se consideran las factorizaciones de Cholesky o LU 
# segun el tipo de matriz, factorizacion QR y resolucion por inversa de la 
# matriz A; en todos los casos se calculan los tiempos y la precision de la
# solucion medida como la desviacion estandar del residuo.

# los resultados son almacenados en un archivo .txt

# Corrida para matrices definidas positivas y simetricas
file = fopen('Proyecto corridas1.txt', 'w')
for i=[2:10,20:10:100,200:100:900]
# generar matriz
n=i;
disp('--------------------------------------')
n
A= rand(n)*100;
[Q,R]=qr(A);
v= rand(1,n)*100;
A=Q'*diag(v)*Q;
A=(A'+A)*(1/2);
b= rand(n,1)*100;

# numero de condicion de A
norma= norm(A);
condicio=norma*norm(inv(A))

# validacion de definida positiva
def=7;
for i=1:n
  detA=det(A(1:i,1:i));
  if (detA<0)
    def=1;
    break;
  end
end 

# resoluciones del sistema Ax=b
# Factorizacion cholesky o LU segun aplique
tic;
if(def==1)
  [L,U]=lu(A);
  x1=U\b;
  x2=L\x1;
  disp('Factorizacion LU');

else
  [L]=Cholesky(A);
  x1=L'\b;
  x2=L\x1;
  disp('Factorizacion Cholesky');
end
t1=toc;

# factorizacion QR
tic;
[Q,R]=qr(A);
x3=R\(Q'*b);

t2=toc;

# solucion por inversa
tic;
x=inv(A)*b;
t3=toc;

#precision de los resultados

p1=std((b-A*x2));
p2=std((b-A*x3));
p3=std((b-A*x));

# Almacenamiento de los resultados
matriz= {'Metodo','Tiempo','Precision';'Cholesky', t1,p1;'Fact QR',t2,p2;'Inversa',t3,p3};

fprintf(file, 'N= %d \n', n);
fprintf(file, 'Metodo \t \t Tiempo \t Precision \n', matriz{1,:});
for j=2:4
  fprintf(file, '%s \t %d \t %d \n', matriz{j,:});
end
fprintf(file, '\n cond(A)= %f \n', condicio);
fprintf(file, '--------------------------------------- \n ');

end
fclose(file);


# corrida para matrices aleatorias
file = fopen('Proyecto corridas sim.txt', 'w')
for i=[2:10,20:10:100,200:100:900]
# generar matriz
n=i;
disp('--------------------------------------')
n
A= rand(n)*100;
b= rand(n,1)*100;

# numero de condicion de la matriz A
norma= norm(A);
condicio=norma*norm(inv(A))

# validacion de matriz definida positiva
def=7;
for i=1:n
  detA=det(A(1:i,1:i));
  if (detA<0)
    def=1;
    break;
  end
end 

# resolucion del sistema Ax=b
# Factorizacion cholesky o Lu segun aplique
tic;
if(def==1)
  [L,U]=lu(A);
  x1=U\b;
  x2=L\x1;
  disp('Factorizacion LU');
dis='Fact LU';
else
  [L]=Cholesky(A);
  x1=L'\b;
  x2=L\x1;
  disp('Factorizacion Cholesky');
end
t1=toc;

# factorizacion QR
tic;
[Q,R]=qr(A);
x3=R\(Q'*b);
t2=toc;

# solucion por medio de inversa
tic;
x=inv(A)*b;
t3=toc;

# precision de los resultados

p1=std((b-A*x2));
p2=std((b-A*x3));
p3=std((b-A*x));

# almacenamiento de datos
matriz= {'Metodo','Tiempo','Precision';dis, t1,p1;'Fact QR',t2,p2;'Inversa',t3,p3};

fprintf(file, 'N= %d \n', n);
fprintf(file, 'Metodo \t \t Tiempo \t Precision \n', matriz{1,:});
for j=2:4
  fprintf(file, '%s \t %d \t %d \n', matriz{j,:});
end
fprintf(file, '\n cond(A)= %f \n', condicio);
fprintf(file, '--------------------------------------- \n ');

end
fclose(file);

