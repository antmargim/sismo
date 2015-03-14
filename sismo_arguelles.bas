cls
ndP=2        'Numero de porticos en la misma direccion
qp=510       'kg/m2 carga permnente
qu=300       'carga de uso
cp.qp=1      'coeficiente de ponderacion de la carga permanente
cp.qu=0.6    'coeficienbte de ponderacion de la craga de uso
lon.for=64.8 'longitud del forjado
anc.for=21.6 'anchura del forjado
ndp=8        'numero de pisos
aPb=4        'altura de planta baja
aPi=3        'altura de piso
ab=0.11      'aceleracion sismica basica
p=1          'coeficiente de riesgo
beta=0.27    'coeficiente de respuesta

qTp=(qp*cp.qp +qu*cp.qu)*lon.for*anc.for
qTp=970000
print "Analisis sismico"
print
print "1.- Estimacion de cargas"
print "(";t$(qp,0);"*";d$(cp.qp,1);"*";t$(qu,0);"*";d$(cp.qu,1);")*";t$(lon.for,2);"*";t$(anc.for,2);" = ";v$(qTp,0)

qTc=800000
print
print "2.- Espectro de respuesta"
k=1   'coeficiente de contribucion
C=1.8 'coeficiente de suelo
alfaTo=(3*C-3.8)*(k-1.25)+2.3
print "Ordenada espectral del periodo p.To = ";d$(alfaTo,2)
p.To=0.125*C+0.2*k-0.175
print "Periodo  p.To                       = ";d$(p.To,3);" sg"
p.T1=0.215*k*(5*C-1)/alfaTo
print "Periodo  p.T1                       = ";d$(p.T1,3);" sg"

print
print "3.- Periodos de los modos de vibracion"

for n=1 To 3
 p.T(n)=0.11*ndp/(n*2-1)
 alfaT1(i)=alfaTo*p.T(n)/n
 print "Periodo del modo ";n;" p.T =";d$(p.T(n),3);" sg  "
next n
print

for i=1 to 3    'ordenadas espectrales para el periodo de cada modo
 if p.T(i)<p.To then alfaT(i)=1+(alfaTo-1)*p.T(i)/p.To
 if p.T(i)>=p.To and p.T(i)<=p.T1 then alfaT(i)=alfaTo
 if p.T(i)>p.T1 then alfaT(i)=alfaTo*p.T(i)/p.To
next i
print
print "4.- Factores de distribucion"

pi=4*atn(1)
H.t=aPb
for n=2 to ndp
 H.t=H.t+aPi
next n

print "Altura total = ";H.t
print
print "F = Fuerzas estaticas    V = Cortante por planta"
print "h.piso Carga c.Forma coefi f.distri  c.sismico   F   V"
print "-------------------------------------------------------"

for i=1 to 3
print "            Modo ";i
coef1=0
coef2=0
for n=1 to ndp
 h.piso(n)=aPi*(n-1)+aPb
 qp(n)=qTp:if n=ndp then qp(n)=qTc
 co.for(n)=sin((2*i-1)*pi*h.piso(n)/2/H.t)
 coef1=coef1+qp(n)*co.for(n)
 coef2=coef2+qp(n)*co.for(n)^2
next n

coefi=coef1/coef2
for n=1 to ndp
 nu(n)=coefi*co.for(n)
 s(n)=ab*p*alfaT(i)*beta*nu(n)
 F(n)=s(n)*qp(n)/1000
next n

V(ndp+1,i)=0  'cortantes por planta
for n=ndp to 1 step -1
 V(n,i)=F(n)+V(n+1,i)
next n

for n=1 to ndp
 pepe$="       ":if n=1 then pepe$=d$(coefi,4)
 print "  ";d$(h.piso(n),0);" ";t$(qp(n)/1000,0);" ";d$(co.for(n),4);" ";pepe$;" ";d$(nu(n),4);_
     "    ";d$(s(n),4);" ";t$(F(n),2);" ";t$(V(n,i),2)
next n
print ""
next i

for n=1 to ndp   'Cortantes combinados
V.C2(n)=sqr(V(n,1)^2+V(n,2)^2)
V.C3(n)=sqr(V(n,1)^2+V(n,2)^2+V(n,3)^2)
next n
print
for n=1 to ndp   'Fuerzas equivalentes
F.e1(n)=V(n,1)-V(n+1,1)
F.e2(n)=V.C2(n)-V.C2(n+1)
F.e3(n)=V.C3(n)-V.C3(n+1)
next n
print
for n=1 to ndp  'Fuerzas por portico
F1(n)=F.e1(n)/ndP
F2(n)=F.e2(n)/ndP
F3(n)=F.e3(n)/ndP
next n
print
print "V.C=Cortante combinado  F=Fuerza por portico  F.e=Fuerza equivalentes"
print "            1 modo                2 modos                3 modos"
print " n     V.C   F.e     F       V.C   F.e     F         V.C   F.e     F"
print "-----------------------------------------------------------------------"
for n=1 to ndp
print d$(n,0);"  ";t$(V(n,1),2);" ";t$(F.e1(n),2);" ";t$(F1(n),2);"   ";_
                  t$(V.C2(n),2);" ";t$(F.e2(n),2);" ";t$(F2(n),2);"   ";_
                  t$(V.C3(n),2);" ";t$(F.e3(n),2);" ";t$(F3(n),2)
next n
wait


'...........................................................
sub de n0,k0
 global en,Le$
 D$=""
 E$=""
 for n=1 to n0:D$=D$+"#":next
 For n=1 to en:E$=E$+"#":next
 Le$=using (E$+"."+D$,k0)
end sub
'.................
function u$(kk,n1)
en=1:call de n1,kk:u$=Le$
end function

function d$(kk,n2)
en=2:call de n2,kk:d$=Le$
end function

function t$(kk,n3)
en=3:call de n3,kk:t$=Le$
end function

function q$(kk,n4)
en=4:call de n4,kk:q$=Le$
end function

function f$(kk,n5)
en=5:call de n5,kk:f$=Le$
end function

function s$(kk,n6)
en=6:call de n6,kk:s$=Le$
end function

function v$(kk,n7)
en=7:call de n7,kk:v$=Le$
end function

function o$(kk,n8)
en=8:call de n8,kk:o$=Le$
end function
