import numpy as np
from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *
from operator import add

def jednacine(original,slika,T1,T2):
	M=[]
	M.append(np.subtract(np.dot(original[1],T1[2]),np.dot(original[2],T1[1])))
	M.append(np.add(np.dot(-1*original[0],T1[2]),np.dot(original[2],T1[0])))
	M.append(np.subtract(np.dot(slika[1],T2[2]),np.dot(slika[2],T2[1])))
	M.append(np.add(np.dot(-1*slika[0],T2[2]),np.dot(slika[2],T2[0])))
	return M

def Matrica(v):
	M=[]
	v1=[0,-1*v[2],v[1]]
	v2=[v[2],0,-1*v[0]]
	v3=[-1*v[1],v[0],0]
	M.append(v1)
	M.append(v2)
	M.append(v3)
	return M

def Red(a,b):
	R=[a[0]*b[0],a[1]*b[0],a[2]*b[0],a[0]*b[1],a[1]*b[1],a[2]*b[1],a[0]*b[2],a[1]*b[2],a[2]*b[2]]
	return R

def DLT(Originali,Slike,n):
	P=[]
	for i in range(0,n):
		R=Red(Originali[i],Slike[i])
		P.append(R)

	P1=np.linalg.svd(P)
	Zadnja=P1[2]
	Trazena=Zadnja[8]
	Matrica=[]

	Prvi=[]
	Prvi.append(Trazena[0])
	Prvi.append(Trazena[1])
	Prvi.append(Trazena[2])

	Drugi=[]
	Drugi.append(Trazena[3])
	Drugi.append(Trazena[4])
	Drugi.append(Trazena[5])

	Treci=[]
	Treci.append(Trazena[6])
	Treci.append(Trazena[7])
	Treci.append(Trazena[8])

	Matrica.append(Prvi)
	Matrica.append(Drugi)
	Matrica.append(Treci)

	return(Matrica)

def izracunaj():

	Leva1=[332.0,72.0,1.0]
	Leva2=[502.0,56.0,1.0]
	Leva3=[718.0,167.0,1.0]
	Leva4=[538.0,188.0,1.0]
	Leva5=[329.0,292.0,1.0]
	Leva7=[714.0,400.0,1.0]
	Leva8=[536.0,428.0,1.0]
	Leva9=[263.0,339.0,1.0]
	Leva11=[775.0,365.0,1.0]
	Leva12=[313.0,416.0,1.0]
	Leva13=[266.0,587.0,1.0]
	Leva15=[768.0,616.0,1.0]
	Leva16=[315.0,668.0,1.0]
	Leva17=[90.0,630.0,1.0]
	Leva19=[926.0,599.0,1.0]
	Leva20=[700.0,780.0,1.0]
	Leva21=[95.0,830.0,1.0]
	Leva23=[920.0,786.0,1.0]
	Leva24=[700.0,990.0,1.0]

	Desna1=[392.0,76.0,1.0]
	Desna2=[563.0,76.0,1.0]
	Desna3=[567.0,196.0,1.0]
	Desna4=[370.0,197.0,1.0]
	Desna7=[566.0,423.0,1.0]
	Desna8=[375.0,421.0,1.0]
	Desna9=[280.0,312.0,1.0]
	Desna10=[715.0,330.0,1.0]
	Desna11=[687.0,403.0,1.0]
	Desna12=[237.0,380.0,1.0]
	Desna14=[713.0,569.0,1.0]
	Desna15=[685.0,641.0,1.0]
	Desna16=[247.0,615.0,1.0]
	Desna17=[122.0,550.0,1.0]
	Desna19=[861.0,654.0,1.0]
	Desna20=[458.0,777.0,1.0]
	Desna21=[128.0,719.0,1.0]
	Desna23=[857.0,839.0,1.0]
	Desna24=[462.0,974.0,1.0]

	Leva=[]
	Desna=[]

	Leva.append(Leva1)
	Leva.append(Leva2)
	Leva.append(Leva3)
	Leva.append(Leva4)
	Leva.append(Leva11)
	Leva.append(Leva12)
	Leva.append(Leva15)
	Leva.append(Leva16)
	Desna.append(Desna1)
	Desna.append(Desna2)
	Desna.append(Desna3)
	Desna.append(Desna4)
	Desna.append(Desna11)
	Desna.append(Desna12)
	Desna.append(Desna15)
	Desna.append(Desna16)
	
	n=len(Leva)
	greska=0.00001

	#matrica dobijena malo drugacijim DLT algoritmom 

	F=DLT(Leva,Desna,n)

	#odredjivanje epipola1

	ind=True
	for i in range(0,n):
		it=np.dot(Desna[i],np.dot(F,np.matrix.transpose(np.asarray(Leva[i]))))
		if(it>greska or it<-greska):
			ind=False
	d=np.linalg.det(F)
	if(d>greska or d<-greska):
		ind=False
	FF=np.linalg.svd(F)
	e1=FF[2][2]
	e1=[e1[0]/e1[2],e1[1]/e1[2],e1[2]/e1[2]]

	#matrica F transponovana

	FT=np.matrix.transpose(np.asarray(F))

	#odredjivanje epipola2

	ind=True
	for i in range(0,n):
		it=np.dot(Leva[i],np.dot(FT,np.matrix.transpose(np.asarray(Desna[i]))))
		if(it>greska or it<-greska):
			ind=False
	d1=np.linalg.det(FT)
	if(d1>greska or d1<-greska):
		ind=False
	FFT=np.linalg.svd(FT)
	e2=FFT[2][2]
	e2=[e2[0]/e2[2],e2[1]/e2[2],e2[2]/e2[2]]

	#odredjivanje skrivenih tacaka
	
	#mala	

	MaliNedogled1=np.cross(np.cross(Leva4,Leva1),np.cross(Leva3,Leva2))
	MaliNedogled2=np.cross(np.cross(Leva4,Leva3),np.cross(Leva8,Leva7))
	
	Leva6=np.cross(np.cross(Leva7,MaliNedogled1),np.cross(Leva5,MaliNedogled2))
	Leva6=[Leva6[0]/Leva6[2],Leva6[1]/Leva6[2],Leva6[2]/Leva6[2]]

	MaliNedogled1d=np.cross(np.cross(Desna4,Desna1),np.cross(Desna3,Desna2))	
	MaliNedogled2d=np.cross(np.cross(Desna4,Desna8),np.cross(Desna3,Desna7))
	MaliNedogled3d=np.cross(np.cross(Desna4,Desna3),np.cross(Desna8,Desna7))	

	Desna5=np.cross(np.cross(MaliNedogled1d,Desna8),np.cross(MaliNedogled2d,Desna1))
	Desna5=[Desna5[0]/Desna5[2],Desna5[1]/Desna5[2],Desna5[2]/Desna5[2]]
	Desna6=np.cross(np.cross(Desna7,MaliNedogled1d),np.cross(Desna5,MaliNedogled3d))
	Desna6=[Desna6[0]/Desna6[2],Desna6[1]/Desna6[2],Desna6[2]/Desna6[2]]

	#srednja
	
	SrednjiNedogled1=np.cross(np.cross(Leva12,Leva11),np.cross(Leva16,Leva15))
	SrednjiNedogled2=np.cross(np.cross(Leva12,Leva9),np.cross(Leva16,Leva13))
	SrednjiNedogled3=np.cross(np.cross(Leva12,Leva16),np.cross(Leva11,Leva15))	
	SrednjiNedogled1d=np.cross(np.cross(Desna16,Desna15),np.cross(Desna9,Desna10))
	SrednjiNedogled2d=np.cross(np.cross(Desna14,Desna15),np.cross(Desna10,Desna11))
	SrednjiNedogled3d=np.cross(np.cross(Desna12,Desna16),np.cross(Desna10,Desna14))

	Leva10=np.cross(np.cross(SrednjiNedogled1,Leva9),np.cross(SrednjiNedogled2,Leva11))
	Leva10=[Leva10[0]/Leva10[2],Leva10[1]/Leva10[2],Leva10[2]/Leva10[2]]
	Leva14=np.cross(np.cross(SrednjiNedogled1,Leva13),np.cross(Leva15,SrednjiNedogled2))
	Leva14=[Leva14[0]/Leva14[2],Leva14[1]/Leva14[2],Leva14[2]/Leva14[2]]
	Desna13=np.cross(np.cross(SrednjiNedogled1d,Desna14),np.cross(SrednjiNedogled2d,Desna16))
	Desna13=[Desna13[0]/Desna13[2],Desna13[1]/Desna13[2],Desna13[2]/Desna13[2]]

	#velika
	
	VelikiNedogled1=np.cross(np.cross(Leva17,Leva20),np.cross(Leva21,Leva24))
	VelikiNedogled2=np.cross(np.cross(Leva20,Leva19),np.cross(Leva24,Leva23))

	Leva18=np.cross(np.cross(VelikiNedogled2,Leva17),np.cross(VelikiNedogled1,Leva19))
	Leva18=[Leva18[0]/Leva18[2],Leva18[1]/Leva18[2],Leva18[2]/Leva18[2]]
	Leva22=np.cross(np.cross(VelikiNedogled2,Leva21),np.cross(VelikiNedogled1,Leva23))
	Leva22=[Leva22[0]/Leva22[2],Leva22[1]/Leva22[2],Leva22[2]/Leva22[2]]

	VelikiNedogled1d=np.cross(np.cross(Desna17,Desna20),np.cross(Desna21,Desna24))
	VelikiNedogled2d=np.cross(np.cross(Desna20,Desna19),np.cross(Desna24,Desna23))

	Desna18=np.cross(np.cross(VelikiNedogled1d,Desna19),np.cross(VelikiNedogled2d,Desna17))
	Desna18=[Desna18[0]/Desna18[2],Desna18[1]/Desna18[2],Desna18[2]/Desna18[2]]
	Desna22=np.cross(np.cross(VelikiNedogled1d,Desna23),np.cross(VelikiNedogled2d,Desna21))
	Desna22=[Desna22[0]/Desna22[2],Desna22[1]/Desna22[2],Desna22[2]/Desna22[2]]
	
	#triangulacija

	T1=[
	[1,0,0,0],
	[0,1,0,0],
	[0,0,1,0]]

	E2=Matrica(e2)

	T2=[]
	temp=np.matrix.transpose(np.matmul(E2,F))
	T2.append(temp[0])
	T2.append(temp[1])
	T2.append(temp[2])
	T2.append(e2)
	T2=np.matrix.transpose(np.asarray(T2))

	NovaLeva=[]
	NovaDesna=[]	

	NovaLeva.append(Leva1)
	NovaLeva.append(Leva2)
	NovaLeva.append(Leva3)
	NovaLeva.append(Leva4)
	NovaLeva.append(Leva5)
	NovaLeva.append(Leva6)
	NovaLeva.append(Leva7)
	NovaLeva.append(Leva8)
	NovaLeva.append(Leva9)
	NovaLeva.append(Leva10)
	NovaLeva.append(Leva11)
	NovaLeva.append(Leva12)
	NovaLeva.append(Leva13)
	NovaLeva.append(Leva14)
	NovaLeva.append(Leva15)	
	NovaLeva.append(Leva16)
	NovaLeva.append(Leva17)
	NovaLeva.append(Leva18)
	NovaLeva.append(Leva19)
	NovaLeva.append(Leva20)
	NovaLeva.append(Leva21)
	NovaLeva.append(Leva22)
	NovaLeva.append(Leva23)
	NovaLeva.append(Leva24)
	
	NovaDesna.append(Desna1)
	NovaDesna.append(Desna2)
	NovaDesna.append(Desna3)
	NovaDesna.append(Desna4)
	NovaDesna.append(Desna5)
	NovaDesna.append(Desna6)
	NovaDesna.append(Desna7)
	NovaDesna.append(Desna8)
	NovaDesna.append(Desna9)
	NovaDesna.append(Desna10)
	NovaDesna.append(Desna11)
	NovaDesna.append(Desna12)
	NovaDesna.append(Desna13)
	NovaDesna.append(Desna14)
	NovaDesna.append(Desna15)
	NovaDesna.append(Desna16)
	NovaDesna.append(Desna17)
	NovaDesna.append(Desna18)
	NovaDesna.append(Desna19)
	NovaDesna.append(Desna20)
	NovaDesna.append(Desna21)
	NovaDesna.append(Desna22)
	NovaDesna.append(Desna23)
	NovaDesna.append(Desna24)	

	#rekonstrukcija

	n=len(NovaLeva)
	rekonstruisaneTacke=[]
	for i in range(0,n):
		M=[]
		M=jednacine(NovaLeva[i],NovaDesna[i],T1,T2)
		SVD=np.linalg.svd(M)
		tacka=np.matrix.transpose(np.asarray(SVD[2][3]))
		rekonstruisaneTacke.append(tacka)

	#konacne koordinate sa dovoljno skaliranom z koordinatom

	konacneTacke=[]
	for i in range(0,n):
		konacneTacke.append([rekonstruisaneTacke[i][0]/rekonstruisaneTacke[i][3]*10,rekonstruisaneTacke[i][1]/rekonstruisaneTacke[i][3]*10,rekonstruisaneTacke[i][2]/rekonstruisaneTacke[i][3]*(400)*10])
	for i in range(0,len(konacneTacke)):
		print(NovaLeva[i])
		print(NovaDesna[i])
	for i in range(0,len(konacneTacke)):
		print("\n")
		print(i+1,konacneTacke[i])
	return konacneTacke

def iscrtaj(konacneTacke):

#mala	
#	glBegin(GL_LINES)
#	glColor3f(1,0,0)	
#	glVertex3f(konacneTacke[0][0],konacneTacke[0][1],konacneTacke[0][2])
#	glVertex3f(konacneTacke[1][0],konacneTacke[1][1],konacneTacke[1][2])
#	glEnd()
#	glBegin(GL_LINES)
#	glColor3f(1,0,0)
#	glVertex3f(konacneTacke[1][0],konacneTacke[1][1],konacneTacke[1][2])
#	glVertex3f(konacneTacke[2][0],konacneTacke[2][1],konacneTacke[2][2])
#	glEnd()
	glBegin(GL_LINES)
	glColor3f(1,0,0)
	glVertex3f(konacneTacke[2][0],konacneTacke[2][1],konacneTacke[2][2])
	glVertex3f(konacneTacke[3][0],konacneTacke[3][1],konacneTacke[3][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(1,0,0)	
	glVertex3f(konacneTacke[3][0],konacneTacke[3][1],konacneTacke[3][2])
	glVertex3f(konacneTacke[0][0],konacneTacke[0][1],konacneTacke[0][2])
	glEnd()	

	glBegin(GL_LINES)
	glColor3f(1,0,0)
	glVertex3f(konacneTacke[4][0],konacneTacke[4][1],konacneTacke[4][2])
	glVertex3f(konacneTacke[5][0],konacneTacke[5][1],konacneTacke[5][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(1,0,0)	
	glVertex3f(konacneTacke[5][0],konacneTacke[5][1],konacneTacke[5][2])
	glVertex3f(konacneTacke[6][0],konacneTacke[6][1],konacneTacke[6][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(1,0,0)
	glVertex3f(konacneTacke[6][0],konacneTacke[6][1],konacneTacke[6][2])
	glVertex3f(konacneTacke[7][0],konacneTacke[7][1],konacneTacke[7][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(1,0,0)
	glVertex3f(konacneTacke[7][0],konacneTacke[7][1],konacneTacke[7][2])
	glVertex3f(konacneTacke[4][0],konacneTacke[4][1],konacneTacke[4][2])
	glEnd()

	glBegin(GL_LINES)
	glColor3f(1,0,0)		
	glVertex3f(konacneTacke[0][0],konacneTacke[0][1],konacneTacke[0][2])
	glVertex3f(konacneTacke[4][0],konacneTacke[4][1],konacneTacke[4][2])
	glEnd()	
	glBegin(GL_LINES)
	glColor3f(1,0,0)		
	glVertex3f(konacneTacke[2][0],konacneTacke[2][1],konacneTacke[2][2])
	glVertex3f(konacneTacke[6][0],konacneTacke[6][1],konacneTacke[6][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(1,0,0)		
	glVertex3f(konacneTacke[3][0],konacneTacke[3][1],konacneTacke[3][2])
	glVertex3f(konacneTacke[7][0],konacneTacke[7][1],konacneTacke[7][2])
	glEnd()
#	glBegin(GL_LINES)
#	glColor3f(1,0,0)		
#	glVertex3f(konacneTacke[1][0],konacneTacke[1][1],konacneTacke[1][2])
#	glVertex3f(konacneTacke[5][0],konacneTacke[5][1],konacneTacke[5][2])
#	glEnd()

#srednja
	glBegin(GL_LINES)
	glColor3f(0,1,0)	
	glVertex3f(konacneTacke[8][0],konacneTacke[8][1],konacneTacke[8][2])
	glVertex3f(konacneTacke[9][0],konacneTacke[9][1],konacneTacke[9][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(konacneTacke[9][0],konacneTacke[9][1],konacneTacke[9][2])
	glVertex3f(konacneTacke[10][0],konacneTacke[10][1],konacneTacke[10][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(konacneTacke[10][0],konacneTacke[10][1],konacneTacke[10][2])
	glVertex3f(konacneTacke[11][0],konacneTacke[11][1],konacneTacke[11][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)	
	glVertex3f(konacneTacke[11][0],konacneTacke[11][1],konacneTacke[11][2])
	glVertex3f(konacneTacke[8][0],konacneTacke[8][1],konacneTacke[8][2])
	glEnd()

	glBegin(GL_LINES)
	glColor3f(0,1,0)	
	glVertex3f(konacneTacke[12][0],konacneTacke[12][1],konacneTacke[12][2])
	glVertex3f(konacneTacke[13][0],konacneTacke[13][1],konacneTacke[13][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(konacneTacke[13][0],konacneTacke[13][1],konacneTacke[13][2])
	glVertex3f(konacneTacke[14][0],konacneTacke[14][1],konacneTacke[14][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(konacneTacke[14][0],konacneTacke[14][1],konacneTacke[14][2])
	glVertex3f(konacneTacke[15][0],konacneTacke[15][1],konacneTacke[15][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)	
	glVertex3f(konacneTacke[15][0],konacneTacke[15][1],konacneTacke[15][2])
	glVertex3f(konacneTacke[12][0],konacneTacke[12][1],konacneTacke[12][2])
	glEnd()

	glBegin(GL_LINES)
	glColor3f(0,1,0)	
	glVertex3f(konacneTacke[8][0],konacneTacke[8][1],konacneTacke[8][2])
	glVertex3f(konacneTacke[12][0],konacneTacke[12][1],konacneTacke[12][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(konacneTacke[11][0],konacneTacke[11][1],konacneTacke[11][2])
	glVertex3f(konacneTacke[15][0],konacneTacke[15][1],konacneTacke[15][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(konacneTacke[10][0],konacneTacke[10][1],konacneTacke[10][2])
	glVertex3f(konacneTacke[14][0],konacneTacke[14][1],konacneTacke[14][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)	
	glVertex3f(konacneTacke[9][0],konacneTacke[9][1],konacneTacke[9][2])
	glVertex3f(konacneTacke[13][0],konacneTacke[13][1],konacneTacke[13][2])
	glEnd()	

#velika
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[16][0],konacneTacke[16][1],konacneTacke[16][2])
	glVertex3f(konacneTacke[17][0],konacneTacke[17][1],konacneTacke[17][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[17][0],konacneTacke[17][1],konacneTacke[17][2])
	glVertex3f(konacneTacke[18][0],konacneTacke[18][1],konacneTacke[18][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[18][0],konacneTacke[18][1],konacneTacke[18][2])
	glVertex3f(konacneTacke[19][0],konacneTacke[19][1],konacneTacke[19][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[19][0],konacneTacke[19][1],konacneTacke[19][2])
	glVertex3f(konacneTacke[16][0],konacneTacke[16][1],konacneTacke[16][2])
	glEnd()

	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[20][0],konacneTacke[20][1],konacneTacke[20][2])
	glVertex3f(konacneTacke[21][0],konacneTacke[21][1],konacneTacke[21][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[21][0],konacneTacke[21][1],konacneTacke[21][2])
	glVertex3f(konacneTacke[22][0],konacneTacke[22][1],konacneTacke[22][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[22][0],konacneTacke[22][1],konacneTacke[22][2])
	glVertex3f(konacneTacke[23][0],konacneTacke[23][1],konacneTacke[23][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[23][0],konacneTacke[23][1],konacneTacke[23][2])
	glVertex3f(konacneTacke[20][0],konacneTacke[20][1],konacneTacke[20][2])
	glEnd()

	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[16][0],konacneTacke[16][1],konacneTacke[16][2])
	glVertex3f(konacneTacke[20][0],konacneTacke[20][1],konacneTacke[20][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[17][0],konacneTacke[17][1],konacneTacke[17][2])
	glVertex3f(konacneTacke[21][0],konacneTacke[21][1],konacneTacke[21][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[18][0],konacneTacke[18][1],konacneTacke[18][2])
	glVertex3f(konacneTacke[22][0],konacneTacke[22][1],konacneTacke[22][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[19][0],konacneTacke[19][1],konacneTacke[19][2])
	glVertex3f(konacneTacke[23][0],konacneTacke[23][1],konacneTacke[23][2])
	glEnd()

def on_display():
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	konacneTacke=izracunaj()
	iscrtaj(konacneTacke)
	glutSwapBuffers()

def main():
	glutInit()
	glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH|GLUT_DOUBLE)
	glutInitWindowSize(800, 800)
	glutCreateWindow("3Drekonstrukcija")
	glViewport(0, 0, 800, 800)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	gluPerspective(60, 800/800, 1, 10000 )
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	gluLookAt(300, 150, 300, 0, 0, 0, 0, 1, 0)
	glClearColor(0, 0, 0, 0)
	glEnable(GL_DEPTH_TEST)
	glutDisplayFunc(on_display)
	glutMainLoop()

if __name__ == "__main__":
	main()
