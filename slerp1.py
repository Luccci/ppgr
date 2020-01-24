import numpy as np
import math
from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

c1=[3,4,5]
fiP=math.pi
tetaP=math.pi/4
psiP=math.pi/3
c2=[-4,5,3]
fiK=math.pi/6
tetaK=math.pi/1
psiK=math.pi/2
staticni=[]
pocetak=[]
kraj=[]

t=0
tm=1
anim=False
mili=10
pomeraj=0.01

def Euler2A(fi,teta,psi):
	Rx=[[1,0,0],[0,math.cos(fi),math.sin(fi)*(-1)],[0,math.sin(fi),math.cos(fi)]]
	Ry=[[math.cos(teta),0,math.sin(teta)],[0,1,0],[math.sin(teta)*(-1),0,math.cos(teta)]]
	Rz=[[math.cos(psi),math.sin(psi)*(-1),0],[math.sin(psi),math.cos(psi),0],[0,0,1]]
	A=np.linalg.multi_dot([Rz,Ry,Rx])
	return A

def AxisAngle(A):
	nulaZ=0.005
	jedanG=1.005
	jedanD=0.995
	if((A[0][1]>-nulaZ and A[0][1]<nulaZ)and(A[0][2]>-nulaZ and A[0][2]<nulaZ)
        and(A[1][0]>-nulaZ and A[1][0]<nulaZ)and(A[1][2]>-nulaZ and A[1][2]<nulaZ)
        and(A[2][0]>-nulaZ and A[2][0]<nulaZ)and(A[2][1]>-nulaZ and A[2][1]<nulaZ)):
		if(abs(A[0][0]-A[1][1])<nulaZ and abs(A[1][1]-A[2][2]<nulaZ) and abs(A[0][0]-A[2][2]<nulaZ)):
			print("Jedinicna je Matrica A u 2.oj funkciji")
			return
	ort=np.dot(A,np.matrix.transpose(np.asarray(A)))
	ind=0
	if((ort[0][0]<jedanG and ort[0][0]>jedanD)and(ort[1][1]<jedanG and ort[1][1]>jedanD)and(ort[2][2]<jedanG and ort[2][2]>jedanD)
        and(ort[0][1]>-nulaZ and ort[0][1]<nulaZ)and(ort[0][2]>-nulaZ and ort[0][2]<nulaZ)
        and(ort[1][0]>-nulaZ and ort[1][0]<nulaZ)and(ort[1][2]>-nulaZ and ort[1][2]<nulaZ)
        and(ort[2][0]>-nulaZ and ort[2][0]<nulaZ)and(ort[2][1]>-nulaZ and ort[2][1]<nulaZ)):
		ind=0
	else:
		ind=1
	if(ind):
		print("Matrica u 2.oj funkciji nije ortogonalna")
		return
	det=np.linalg.det(A)
	if(det>jedanG or det<jedanD):
		print("determinanta u drugoj funkciji nije 1")
		return
	p=[]
	ind=1
	B=[[A[0][0]-1,A[0][1]+0,A[0][2]+0],[A[1][0]+0,A[1][1]-1,A[1][2]+0],[A[2][0]+0,A[2][1]+0,A[2][2]-1]]
	if((B[0][0]>-nulaZ and B[0][0]<nulaZ)and(B[0][1]>-nulaZ and B[0][1]<nulaZ)and(B[0][2]>-nulaZ and B[0][2]<nulaZ)):
		B[0][0]=B[0][0]+1
		p=np.linalg.solve(B,[1,0,0])
		ind=0
	if((B[1][0]>-nulaZ and B[1][0]<nulaZ)and(B[1][1]>-nulaZ and B[1][1]<nulaZ)and(B[1][2]>-nulaZ and B[1][2]<nulaZ)):
		B[1][1]=B[1][1]+1
		p=np.linalg.solve(B,[0,1,0])
		ind=0
	if((B[2][0]>-nulaZ and B[2][0]<nulaZ)and(B[2][1]>-nulaZ and B[2][1]<nulaZ)and(B[2][2]>-nulaZ and B[2][2]<nulaZ)):
		B[2][2]=B[2][2]+1
		p=np.linalg.solve(B,[0,0,1])
		ind=0
	if(ind):
		B[0][0]=B[0][0]+1
		p=np.linalg.solve(B,[1,0,0])
	i=math.sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
	pj=[p[0]/i,p[1]/i,p[2]/i]
	if(pj[0]==0 and pj[1]==0 and pj[2]==0 ):
		print("vektor p u drugoj funkciji je 0")
		return
	e1=[1,0,0]
	e2=[0,1,0]
	if(pj[0]==0):
		k=e1
	else:
		k=e2
	u=np.cross(pj,k)
	i=math.sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])
	u=[u[0]/i,u[1]/i,u[2]/i]
	u1=np.dot(A,u)
	i=math.sqrt(u1[0]*u1[0]+u1[1]*u1[1]+u1[2]*u1[2])
	u1=[u1[0]/i,u1[1]/i,u1[2]/i]
	fi=math.acos(np.dot(u,u1))
	if(np.linalg.det([pj,u,u1])<0):
		pj[0]=-1*pj[0]
		pj[1]=-1*pj[1]
		pj[2]=-1*pj[2]
	return pj,fi

def Rodrigez(p,fi):
	pIntenzitet=math.sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
	pj=[1.0*p[0]/pIntenzitet,1.0*p[1]/pIntenzitet,1.0*p[2]/pIntenzitet]
	P=np.matrix.transpose(np.asarray(pj))
	PT=np.asarray(pj)
	ppt=np.outer(P,PT)
	E=[[1,0,0],[0,1,0],[0,0,1]]
	Px=np.asarray([[0,pj[2]*(-1),pj[1]],[pj[2],0,pj[0]*(-1)],[pj[1]*(-1),pj[0],0]])
	Rp=np.add(np.add(ppt,math.cos(fi)*(np.subtract(E,ppt))),math.sin(fi)*Px)
	return Rp

def A2Euler(A):
	fi=0
	teta=0
	psi=0
	if(A[2][0]<1):
		if(A[2][0]>-1):
			psi=math.atan2(A[1][0],A[0][0])
			teta=math.asin(-1*(A[2][0]))
			fi=math.atan2(A[2][1],A[2][2])
		else:
			psi=math.atan2(-1*A[0][1],A[1][1])
			teta=math.pi/2
			fi=0
	else:
		psi=math.atan2(-1*A[0][1],A[1][1])
		teta=-1*math.pi/2
		fi=0
	return fi,teta,psi

def AxisAngle2Q(p,fi):
	if(fi==0):
		return 1
	w=math.cos(fi/2)
	i=math.sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
	pj=[p[0]/i,p[1]/i,p[2]/i]
	s=math.sin(fi/2)
	x=pj[0]*s
	y=pj[1]*s
	z=pj[2]*s
	qj=[x,y,z,w]
	return qj

def Q2AxisAngle(q):
	i=math.sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3])
	qj=[q[0]/i,q[1]/i,q[2]/i,q[3]/i]
	if(qj[3]<0):
		qj[0]=-1*qj[0]
		qj[1]=-1*qj[1]
		qj[2]=-1*qj[2]
		qj[3]=-1*qj[3]
	fi=2*math.acos(qj[3])
	if(qj[3]*qj[3]==1):
		p=[1,0,0]
	else:
		s=math.sin(fi/2)
		p=[qj[0]/s,qj[1]/s,qj[2]/s]
	return [p,fi]

def izracunajZaSlerp(fiP,tetaP,psiP,fiK,tetaK,psiK):
	pocetak1=Euler2A(fiP,tetaP,psiP)
	kraj1=Euler2A(fiK,tetaK,psiK)
	pocetak2=AxisAngle(pocetak1)
	kraj2=AxisAngle(kraj1)
	q1=AxisAngle2Q(pocetak2[0],pocetak2[1])
	q2=AxisAngle2Q(kraj2[0],kraj2[1])
	return q1,q2
	
def Slerp(q1,q2,t,tm):
	kosinusO=np.dot(q1,q2)
	if(kosinusO<0):
		q1=np.dot(-1,q1)
		kosinusO=-kosinusO
	if(kosinusO>0.95):
		return q1
	ugaO=math.acos(kosinusO)
	donji=math.sin(ugaO)
	qs=np.dot(math.sin(ugaO*(1-t/tm))/donji,q1)+np.dot(math.sin(ugaO*t/tm)/donji,q2)
	intenzitet=math.sqrt(qs[0]*qs[0]+qs[1]*qs[1]+qs[2]*qs[2]+qs[3]*qs[3])
	qs[0]=qs[0]/intenzitet
	qs[1]=qs[1]/intenzitet
	qs[2]=qs[2]/intenzitet
	qs[3]=qs[3]/intenzitet
	return qs

def on_display():
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	
	global t
	global tm
	global staticni
	global pocetak
	global kraj
	global c1
	global c2

	v=np.subtract(c2,c1)
	v1=[]
	
	v1.append(v[0]*t+c1[0])
	v1.append(v[1]*t+c1[1])
	v1.append(v[2]*t+c1[2])
	trenutniq=Slerp(staticni[0],staticni[1],t,tm)
	trenutni=Q2AxisAngle(trenutniq)

	r=Rodrigez(trenutni[0],trenutni[1])
	a=A2Euler(r)
	print(a)

	koordinatniPocetak()
	kocka(pocetak[0],pocetak[1],c1)	
	kocka(kraj[0],kraj[1],c2)
	kocka(trenutni[0],trenutni[1],v1)	

	glutSwapBuffers()

def on_reshape(w, h):
	w=800
	h=800

def on_keyboard(key, x, y):
	global mili
	global t
	global anim
	if key==b'k':
		if(anim):
			anim=True
		else:
			anim=True
			glutTimerFunc(mili, on_timer ,0)
	elif key==b's':
		anim=False
	elif key==b'r':
		t=0
		glutPostRedisplay()
		
def on_timer(value):
	global mili
	global anim
	global t
	global tm
	global pomeraj
	if(value!=0):
		return
	if(t<tm):
		t=t+pomeraj
	glutPostRedisplay()
	if(anim):
		glutTimerFunc(mili,on_timer,0)

def kocka(p,ugao,v):
	glPushMatrix()
	glTranslate(v[0],v[1],v[2])
	glRotate(ugao*180/math.pi,p[0],p[1],p[2])
	glColor3f(1,1,1)
	glutWireCube(1)
	glBegin(GL_LINES)
	glColor3f(1,0,0)
	glVertex3f(0,0,0)
	glVertex3f(5,0,0)
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(0,0,0)
	glVertex3f(0,5,0)
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(0,0,0)
	glVertex3f(0,0,5)
	glEnd()
	glPopMatrix()

def koordinatniPocetak():
	glBegin(GL_LINES)
	glColor3f(1,0,0)
	glVertex3f(0,0,0)
	glVertex3f(30,0,0)
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,1,0)
	glVertex3f(0,0,0)
	glVertex3f(0,30,0)
	glEnd()
	glBegin(GL_LINES)
	glColor3f(0,0,1)
	glVertex3f(0,0,0)
	glVertex3f(0,0,30)
	glEnd()

def main():
	glutInit()
	glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH|GLUT_DOUBLE)
	glutInitWindowSize(800, 800)
	glutCreateWindow("slerp")
	glViewport(0, 0, 800, 800)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	gluPerspective(60, 800/800, 1, 100 )
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	gluLookAt(5, 20, 15, 0, 0, 0, 0, 1, 0)
	glClearColor(0, 0, 0, 0)
	glEnable(GL_DEPTH_TEST)

	global fiP
	global tetaP
	global psiP
	global fiK
	global tetaK
	global psiK

	global staticni
	staticni=izracunajZaSlerp(fiP,tetaP,psiP,fiK,tetaK,psiK)
	global pocetak
	pocetak=Q2AxisAngle(staticni[0])
	global kraj
	kraj=Q2AxisAngle(staticni[1])

	glutDisplayFunc(on_display)
	glutKeyboardFunc(on_keyboard)
	glutReshapeFunc(on_reshape)
	glutMainLoop()

if __name__ == "__main__":
	main()
