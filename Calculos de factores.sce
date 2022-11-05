

function [Zs]=Calc_de_Factores(n,g,fc,fg, RMGF, RMGG,f,ro,rc,rg)
    //---------------------------------------
    //INPUTS
//n --> num de conductores; 
//g --> num de cables de guardia
//fc--> flecha del conductor; 
//fg-->flecha del cable de guardia; 
//RMGF--> radio medio geometrico de las fases; 
//RMGG-->radio medio geometrico del hilo de guardia;
//f  --> frecuencia de funcionamiento del sistema
//ro --> resistividad del suelo
//rc  --> resistencia del conductor
//rg  --> resistencia del cable de guardia
    //---------------------------------------
    //OUTPUTS
    //MZ-->Mat de impedancias corregidas
    //Z3x3-->Mat de impedancias corregidas y reducidas a 3x3
    //---------------------------------------
    //---------------------------------------
//coordenadas de los conductores
s=n+g
vx=zeros(s)
    for i=2:s
        vx(i)=input("considerando el origen horizontal en el primer conductor de la izq, ingrese las diastancias horizontales a dicho conductor: " )
    end

vy=zeros(s)
    for i=1:s
        vy(i)=input("ingrese la altura de el primer conductor: " )
    end
    
    Mc=[vx,vy]
    disp("la matriz de coordenadas",Mc)
    Mcco=Mc
    //correccion por flecha
    for i=1:s
        if i<=n then
        Mcco(i,2)=Mc(i,2)-0.7*fc
        else
        Mcco(i,2)=Mc(i,2)-0.7*fg
        end
    end
    disp ("La matriz Mcco (Mat de coord corregidas por flecha) es:", Mcco)
//-----calculo los vectores con las coordenadas
Md=Dis(Mcco,s)
    disp("la matriz de vectores",Md)
//-----calculo las distancias entre vectores (que serian las distancias entre cables)
a=size(Md(:,1))
b=a(1)
    for k=1:b
        v=Md(k,:)
        d=norm(v)
        vd(k)=d
    end
    disp("las distancias entre los conductores es: ", vd)
    //la Matriz dp muestra la matriz de los vectores formados por las coordenadas y en la 3er columna las distancias de los mismos
    dp=[Md(:,:) vd]
    //--------------------------------------------------


Dv=DisImgVec2(Mcco,vd,s)
disp("las distancias entre los conductores y la imagen de los vecinos es: ", Dv)
//Dv muestra la matriz de los vectores formados por las coordenadas de las imagenes de los vecinos 

[Mf,X]=MatrizdeFactores(Mcco,Dv,vd,n,g, RMGF, RMGG)
disp("la matriz de factores es: ", Mf)
disp("la matriz de inductancias X es: ", X)

//-------------------------------------------------------
//calculo de los p y titas
[Mtita, Mp]=titasypees(n,g,Md,Mcco,f,ro)
disp("la matriz de factores titas es: ", Mtita)
disp("la matriz de factores p ciquitos es: ", Mp)
//-------------------------------------------------------

//-------------------------------------------------------
//calculo de P y Q
[P,Q]=MPyQ(Mp,Mtita,n,g,f)
disp("la matriz de factores de correccion P es: ", P)
disp("la matriz de factores de correccion Q es: ", Q)
//-------------------------------------------------------
//
//-------------------------------------------------------
//calculo de Matriz de impedancias (MZ)
[MZ]=MatdeImp(P,Q,X,n,g,rc,rg)
disp("La matriz de Impedancias corregidas MZ es:", MZ)
//-------------------------------------------------------

//-------------------------------------------------------
//Reduccion de Matriz de impedancias (MZ) a 3x3 (Z3x3)
[Z3x3]=reducirZ(MZ,n,g)
disp("La matriz de Impedancias de corregidas y reducidas a 3x3 Z3x3 es:", Z3x3)
//-------------------------------------------------------

//-------------------------------------------------------
//CALCULO MATRIZ DE IMPEDANCIA DE SECUENCIAS (Zs)
        //Operador "a"
        oa=-.5+imult(sqrt(3)/2)
        //Matriz de transformacion 
        A=[1 1 1;1 oa^2 oa;1 oa oa^2]
        //Inversa de "A"
        A1=A^-1
Zs=A1*Z3x3*A
disp("la matriz de secuencias Zs es:")
//-------------------------------------------------------



endfunction


//**************************************************************
//**************************************************************
//**************************************************************


//-----------------------funcion q calcula los vectores con la matriz de coordenadas
function Mm=Dis(M,n)
    //M matriz de n filas y 2 columnas, vectores coordenadas; n--> num de filas
    //-----------Cantidad de combinaciones posibles
    comb=factorial(n)/(factorial(2)*factorial(n-2))
    //-------------------------------------

Mm=zeros(comb,2)
t=1
k=1
for i=1:comb
    k=k+1
    for j=1:2
        if k<=n then
            Mm(i,j)=M(t,j)-M(k,j)
        else
            t=t+1
            k=t+1
            Mm(i,j)=M(t,j)-M(k,j)
        end
    end
end
endfunction
//-----------------------------------------------------------



//------------------------------------------------------------
//funcion q calcula los vectores a las imagenes de los vecinos con la matriz de coordenadas EN PROMEDIO
function Dv=DisImgVec2(Mcco,vd,n)
    //Mcco matriz de n filas y 2 columnas, donde solo me interesan las alturas; n--> num de filas
    //vd vector de distancias entre conductores
    //-----------Cantidad de combinaciones posibles
    comb=factorial(n)/(factorial(2)*factorial(n-2))


//-----------------------------------------
    Dv=zeros(comb,1)
    t=1
    k=1
    for i=1:comb
        k=k+1
                if k<=n then
                    Dv(i)=(((((2*Mcco(t,2))^2)+((vd(i))^2))^(1/2))+((((2*Mcco(k,2))^2)+((vd(i))^2))^(1/2)))/2
                else
                    t=t+1
                    k=t+1
                    Dv(i)=(((((2*Mcco(t,2))^2)+((vd(i))^2))^(1/2))+((((2*Mcco(k,2))^2)+((vd(i))^2))^(1/2)))/2
                end 
    end
endfunction
//-----------------------------------------------------------


//-----------------------------------------------------------
//------CALCULO LA MATRIZ DE FACTORES
function [Mf,X]=MatrizdeFactores(M,Dv,v,n,g, RMGF, RMGG)
//n--> num de conduct ; 
//g--> num de hilos de guard; 
//v--> vector distancias entre conductores; 
//Dv-->vector distancias entre img vecinos, 
//M--> matriz con las coordenadas y alturas corregidas; 
//RMGF--> radio medio geometrico de las fases; 
//RMGG-->radio medio geometrico del hilo de guardia;

s=n+g
Mf=zeros(s,s)
for i=1:s
    if i<=n then
        Mf(i,i)=4.6052e-4*log10(2*M(i,2)/RMGF)
    else
        Mf(i,i)=4.6052e-4*log10(2*M(i,2)/RMGG)
        end
end
t=1
for i=1:s
    for j=(i+1):s
        Mf(i,j)=4.6052e-4*log10(Dv(t)/v(t))
        t=t+1
        Mf(j,i)=Mf(i,j)
    end
end
//------------------
//MAtriz X

X=2*50*%pi*Mf
endfunction
//-------------------------------------------------


//-------------------------------------------------
//------CALCULO LA MATRIZ DE TITAS Y P
function [Mtita, Mp]=titasypees(n,g,Mv,Mcco,f,ro)
    //f  --> frecuencia de funcionamiento del sistema
    //ro --> resistividad del suelo
    //n -->num de conductores
    //g-->num de hilos de guardia
    // Mv es una Matriz con las distancias horizontales entre los conductores colocando el origen en el primer conductor de izq a derecha
    // Mcco es una Matriz con las distancias verticales corregidas de los conductores colocando el origen en el suelo

    
//-------------------------------------------------------
//CALCULO LOS p CHIQUITOS
//matriz para rellenar con los valores de p chiquitos
s=n+g
Mp=zeros(s,s)//matriz para rellenar con los valores de p chiquitos
for i=1:s
        Mp(i,i)=5.620e-3*Mcco(i,2)*sqrt(f/ro)
end
t=1
for i=1:s
    for j=(i+1):s
        Mp(i,j)=28.1004e-4*Dv(t)*sqrt(f/ro)
        t=t+1
        Mp(j,i)=Mp(i,j)
    end
end
//-------------------------------------------


//-------------------------------------------
//CALCULO LOS TITAS 
Mtita=zeros(s,s)

t=1
for i=1:s
    for j=(i+1):s
        Mtita(i,j)=atan(abs(Mv(t,1))/(Mcco(i,2)+Mcco(j,2)))
        t=t+1
        Mtita(j,i)=Mtita(i,j)
    end
end
//------------------------------------------


endfunction
//-------------------------------------------------


//-------------------------------------------------
//------CALCULO LA MATRIZ P (MAYUSCULA) Y Q
function [P,Q]=MPyQ(Mp,Mtita,n,g,f)
        //f  --> frecuencia de funcionamiento del sistema
    //n -->num de conductores
    //g-->num de hilos de guardia
    // Mp es una Matriz con los p chiquitos
    // Mtitas es una Matriz con los titas 
    s=n+g
//-------------------------------------------------------
//CALCULO LA MATRIZ DE CORRECCION P (mayuscula)
Pij=zeros(s,s)
a=%pi/8
b=1/(3*sqrt(2))
c=1/16
for i=1:s
    for j=1:s
        Pij(i,j)=a-b*Mp(i,j)*cos(Mtita(i,j))+c*((Mp(i,j))^2)*cos(2*Mtita(i,j))*(.6728-log(2/Mp(i,j)))+c*((Mp(i,j))^2)*sin(2*Mtita(i,j))
    end
end
P=25.134e-4*f*Pij
//---------------------------------------------

//----------------------------------------------
//CALCULO LA MATRIZ DE CORRECCION Q (mayuscula)
Qij=zeros(s,s)
for i=1:s
    for j=1:s
        Qij(i,j)=-.0386+.5*log(2/Mp(i,j))+b*Mp(i,j)*cos(Mtita(i,j))
    end
end
Q=25.134e-4*f*Qij
//----------------------------------------------
endfunction
//--------------------------------------------------------------

//---------------------------------------------------------------
//------CALCULO LA MATRIZ DE IMPEDANCIAS 
function [MZ]=MatdeImp(P,Q,X,n,g,rc,rg)
    //rc  --> resistencia del conductor
    //rg  --> resistencia del cable de guardia
    //n -->num de conductores
    //g-->num de hilos de guardia
    // P Mat de correcion
    // Q Mat de correcion
    s=n+g
//-------------------------------------------------------
//CALCULO LA MATRIZ DE DE RESISTENCIAS R 
R=zeros(s,s)
a=.0036
t0=20
tf=75
for i=1:s
    if i<=n then
        R(i,i)=rc*(1+a*(tf-t0))
    else
        R(i,i)=rg
    end
end
disp("la matriz de resistencias corregidas por Temp es:", R)
//-------------------------------------------------------

//-------------------------------------------------------
ZP=R+P
disp("la matriz de resistencias corregidas con la matriz P es ZP: ", ZP)
ZQ=X+Q
disp("la matriz de inductancias corregidas con la matriz Q es ZQ: ", ZQ)
MZ=ZP+imult(ZQ)
//-------------------------------------------------------
endfunction
//-----------------------------------------------------------

//-----------------------------------------------------------
//------REDUCCION DE LA MATRIZ DE IMPEDANCIAS A 3X3
function [Z3x3]=reducirZ(MZ,n,g)
//n-->num de conductores
//g numde cables de guardia
//MZ-->Matriz Z sin reducir
s=n+g
for i=1:3
    for j=1:3
        Z3x3(i,j)=MZ(i,j)-(MZ(i,s)*MZ(j,s)/MZ(s,s))
    end
end

endfunction
//-----------------------------------------------------------
