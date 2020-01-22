import numpy as np
from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *
from operator import add
r=1000
ugao=0

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
	Leva1=[958.0,38.0,1.0]
	Desna1=[933.0,33.0,1.0]
	Leva2=[1117.0,111.0,1.0]
	Desna2=[1027.0,132.0,1.0]
	Leva3=[874.0,285.0,1.0]
	Desna3=[692.0,223.0,1.0]
	Leva4=[707.0,218.0,1.0]
	Desna4=[595.0,123.0,1.0]
	Leva9=[292.0,569.0,1.0]
	Desna9=[272.0,360.0,1.0]
	Leva10=[770.0,969.0,1.0]
	Desna10=[432.0,814.0,1.0]
	Leva11=[770.0,1465.0,1.0]
	Desna11=[414.0,1284.0,1.0]
	Leva12=[317.0,1057.0,1.0]
	Desna12=[258.0,818.0,1.0]

	Originali=[]
	Slike=[]
	Originali.append(Leva1)
	Originali.append(Leva2)
	Originali.append(Leva3)
	Originali.append(Leva4)
	Originali.append(Leva9)
	Originali.append(Leva10)
	Originali.append(Leva11)
	Originali.append(Leva12)
	Slike.append(Desna1)
	Slike.append(Desna2)
	Slike.append(Desna3)
	Slike.append(Desna4)
	Slike.append(Desna9)
	Slike.append(Desna10)
	Slike.append(Desna11)
	Slike.append(Desna12)

	n=len(Originali)
	greska=0.00001

	#matrica dobijena malo drugacijim DLT algoritmom 

	F=DLT(Originali,Slike,n)
	print(F)

	#odredjivanje epipola1

	ind=True
	for i in range(0,n):
		it=np.dot(Slike[i],np.dot(F,np.matrix.transpose(np.asarray(Originali[i]))))
		if(it>greska or it<-greska):
			ind=False
		else:
			print(it)
	d=np.linalg.det(F)
	print(d)
	if(d>greska or d<-greska):
		ind=False
	print(ind)
	FF=np.linalg.svd(F)
	e1=FF[2][2]
	e1=[e1[0]/e1[2],e1[1]/e1[2],e1[2]/e1[2]]
	print(e1)
	print(FF)
	#matrica F transponovana

	FT=np.matrix.transpose(np.asarray(F))
	print(FT)

	#odredjivanje epipola2

	ind=True
	for i in range(0,n):
		it=np.dot(Originali[i],np.dot(FT,np.matrix.transpose(np.asarray(Slike[i]))))
		if(it>greska or it<-greska):
			ind=False
	d1=np.linalg.det(FT)
	print(d1)
	if(d1>greska or d1<-greska):
		ind=False
	print(ind)
	FFT=np.linalg.svd(FT)
	e2=FFT[2][2]
	e2=[e2[0]/e2[2],e2[1]/e2[2],e2[2]/e2[2]]
	print(e2)
	print(FFT)

	#ostale tacke

	Leva6=[1094.0,536.0,1.0]
	Desna6=[980.0,535.0,1.0]
	Leva7=[862.0,729.0,1.0]
	Desna7=[652.0,638.0,1.0]
	Leva8=[710.0,648.0,1.0]
	Desna8=[567.0,532.0,1.0]
	Leva14=[1487.0,598.0,1.0]
	Desna14=[1303.0,700.0,1.0]
	Leva15=[1462.0,1079.0,1.0]
	Desna15=[1257.0,1165.0,1.0]
	Desna13=[1077.0,269.0,1.0]

	#odredjivanje skrivenih tacaka
	print("OVDE")
	Leva5=np.cross(np.cross(np.cross(np.cross(Leva4,Leva8),np.cross(Leva6,Leva2)),Leva1),
	np.cross(np.cross(np.cross(Leva1,Leva4),np.cross(Leva3,Leva2)),Leva8))
	Leva5=[Leva5[0]/Leva5[2],Leva5[1]/Leva5[2],Leva5[2]/Leva5[2]]
	print(Leva5)

	Desna5=np.cross(np.cross(np.cross(np.cross(Desna4,Desna8),np.cross(Desna6,Desna2)),Desna1),
	np.cross(np.cross(np.cross(Desna1,Desna4),np.cross(Desna3,Desna2)),Desna8))
	Desna5=[Desna5[0]/Desna5[2],Desna5[1]/Desna5[2],Desna5[2]/Desna5[2]]
	print(Desna5)

	Leva13=np.cross(np.cross(Leva14,np.cross(np.cross(Leva10,Leva9),np.cross(Leva11,Leva12))),
	np.cross(Leva9,np.cross(np.cross(Leva10,Leva14),np.cross(Leva11,Leva15))))
	Leva13=[Leva13[0]/Leva13[2],Leva13[1]/Leva13[2],Leva13[2]/Leva13[2]]
	print(Leva13)

	Leva16=np.cross(np.cross(np.cross(np.cross(Leva10,Leva14),np.cross(Leva11,Leva15)),Leva12),
	np.cross(np.cross(np.cross(Leva9,Leva10),np.cross(Leva12,Leva11)),Leva15))
	Leva16=[Leva16[0]/Leva16[2],Leva16[1]/Leva16[2],Leva16[2]/Leva16[2]]
	print(Leva16)

	Desna16=np.cross(np.cross(np.cross(np.cross(Desna10,Desna14),np.cross(Desna11,Desna15)),Desna12),
	np.cross(np.cross(np.cross(Desna9,Desna10),np.cross(Desna12,Desna11)),Desna15))
	Desna16=[Desna16[0]/Desna16[2],Desna16[1]/Desna16[2],Desna16[2]/Desna16[2]]
	print(Desna16)

	#triangulacija

	T1=[
	[1,0,0,0],
	[0,1,0,0],
	[0,0,1,0]]

	E2=Matrica(e2)
	print(E2)

	T2=[]
	temp=np.matrix.transpose(np.matmul(E2,F))
	T2.append(temp[0])
	T2.append(temp[1])
	T2.append(temp[2])
	T2.append(e2)
	T2=np.matrix.transpose(np.asarray(T2))
	print(T2)

	NoviOriginali=[]
	NoveSlike=[]
	NoviOriginali.append(Leva1)
	NoviOriginali.append(Leva2)
	NoviOriginali.append(Leva3)
	NoviOriginali.append(Leva4)
	NoviOriginali.append(Leva5)
	NoviOriginali.append(Leva6)
	NoviOriginali.append(Leva7)
	NoviOriginali.append(Leva8)
	NoviOriginali.append(Leva9)
	NoviOriginali.append(Leva10)
	NoviOriginali.append(Leva11)
	NoviOriginali.append(Leva12)
	NoviOriginali.append(Leva13)
	NoviOriginali.append(Leva14)
	NoviOriginali.append(Leva15)
	NoviOriginali.append(Leva16)
	NoveSlike.append(Desna1)
	NoveSlike.append(Desna2)
	NoveSlike.append(Desna3)
	NoveSlike.append(Desna4)
	NoveSlike.append(Desna5)
	NoveSlike.append(Desna6)
	NoveSlike.append(Desna7)
	NoveSlike.append(Desna8)
	NoveSlike.append(Desna9)
	NoveSlike.append(Desna10)
	NoveSlike.append(Desna11)
	NoveSlike.append(Desna12)
	NoveSlike.append(Desna13)
	NoveSlike.append(Desna14)
	NoveSlike.append(Desna15)
	NoveSlike.append(Desna16)

	#rekonstrukcija

	n=len(NoviOriginali)
	rekonstruisaneTacke=[]
	for i in range(0,n):
		M=[]
		M=jednacine(NoviOriginali[i],NoveSlike[i],T1,T2)
		SVD=np.linalg.svd(M)
		tacka=np.matrix.transpose(np.asarray(SVD[2][3]))
		rekonstruisaneTacke.append(tacka)

	#konacne koordinate sa dovoljno skaliranom z koordinatom

	konacneTacke=[]
	for i in range(0,n):
		konacneTacke.append([rekonstruisaneTacke[i][0]/rekonstruisaneTacke[i][3],rekonstruisaneTacke[i][1]/rekonstruisaneTacke[i][3],rekonstruisaneTacke[i][2]/rekonstruisaneTacke[i][3]*(400)])
	print(rekonstruisaneTacke)
	print(konacneTacke)
	print(NoviOriginali)
	print(NoveSlike)
	return konacneTacke

def iscrtaj(konacneTacke):
#mala
	glBegin(GL_LINES)
	glColor3f(1,0,0)	
	glVertex3f(konacneTacke[0][0],konacneTacke[0][1],konacneTacke[0][2])
	glVertex3f(konacneTacke[1][0],konacneTacke[1][1],konacneTacke[1][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(1,0,0)
	glVertex3f(konacneTacke[1][0],konacneTacke[1][1],konacneTacke[1][2])
	glVertex3f(konacneTacke[2][0],konacneTacke[2][1],konacneTacke[2][2])
	glEnd()
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
	glVertex3f(konacneTacke[1][0],konacneTacke[1][1],konacneTacke[1][2])
	glVertex3f(konacneTacke[5][0],konacneTacke[5][1],konacneTacke[5][2])
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

#velika
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[8][0],konacneTacke[8][1],konacneTacke[8][2])
	glVertex3f(konacneTacke[9][0],konacneTacke[9][1],konacneTacke[9][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[9][0],konacneTacke[9][1],konacneTacke[9][2])
	glVertex3f(konacneTacke[13][0],konacneTacke[13][1],konacneTacke[13][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[13][0],konacneTacke[13][1],konacneTacke[13][2])
	glVertex3f(konacneTacke[12][0],konacneTacke[12][1],konacneTacke[12][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)	
	glVertex3f(konacneTacke[12][0],konacneTacke[12][1],konacneTacke[12][2])
	glVertex3f(konacneTacke[8][0],konacneTacke[8][1],konacneTacke[8][2])
	glEnd()
	
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[11][0],konacneTacke[11][1],konacneTacke[11][2])
	glVertex3f(konacneTacke[10][0],konacneTacke[10][1],konacneTacke[10][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[10][0],konacneTacke[10][1],konacneTacke[10][2])
	glVertex3f(konacneTacke[14][0],konacneTacke[14][1],konacneTacke[14][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[14][0],konacneTacke[14][1],konacneTacke[14][2])
	glVertex3f(konacneTacke[15][0],konacneTacke[15][1],konacneTacke[15][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[15][0],konacneTacke[15][1],konacneTacke[15][2])
	glVertex3f(konacneTacke[11][0],konacneTacke[11][1],konacneTacke[11][2])
	glEnd()
	
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[8][0],konacneTacke[8][1],konacneTacke[8][2])
	glVertex3f(konacneTacke[11][0],konacneTacke[11][1],konacneTacke[11][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[9][0],konacneTacke[9][1],konacneTacke[9][2])
	glVertex3f(konacneTacke[10][0],konacneTacke[10][1],konacneTacke[10][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[14][0],konacneTacke[14][1],konacneTacke[14][2])
	glVertex3f(konacneTacke[13][0],konacneTacke[13][1],konacneTacke[13][2])
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(konacneTacke[15][0],konacneTacke[15][1],konacneTacke[15][2])
	glVertex3f(konacneTacke[12][0],konacneTacke[12][1],konacneTacke[12][2])
	glEnd()


def on_display():
	global r
	global ugao
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	gluLookAt(r*np.cos(ugao), r, np.sin(ugao)*r, 0, 0, 0, 0, 1, 0)
	konacneTacke=izracunaj()
	iscrtaj(konacneTacke)
	glutSwapBuffers()

def on_keyboard(key, x, y):
	global ugao
	global r
	if(key==b'a'):
		ugao=ugao-0.1
		glutPostRedisplay()
	if(key==b'd'):
		ugao=ugao+0.1
		glutPostRedisplay()
	if(key==b'w'):
		r=r-100
		glutPostRedisplay()
	if(key==b's'):
		r=r+100
		glutPostRedisplay()

def main():
	glutInit()
	glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH|GLUT_DOUBLE)
	glutInitWindowSize(800, 800)
	glutCreateWindow("3Drekonstrukcija")
	glViewport(0, 0, 800, 800)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	gluPerspective(60, 800/800, 1, 1000000 )
	glClearColor(0, 0, 0, 0)
	glEnable(GL_DEPTH_TEST)
	glutDisplayFunc(on_display)
	glutKeyboardFunc(on_keyboard)
	glutMainLoop()

if __name__ == "__main__":
	main()
