import matplotlib.pyplot as plt
import numpy as np
s=[]
kin=[]
rho_0=1
rho_e=rho_0
c_Ae = 4
c_A=3


def saus(c,kza):
    a = rho_e*c_Ae**2*np.tanh(kza)+rho_0*c_A**2
    b = rho_e*np.tanh(kza)+rho_0
    return c**2-a/b
def kink(c,kza):
    a=rho_e*c_Ae**2+rho_0*c_A**2*np.tanh(kza)
    b= rho_e+rho_0*np.tanh(kza)
    return c**2-a/b

def root_stepper(f,x,y,dx,w):
    x1 = x ; f1 = f(x,w)
    x2 = x +dx; f2 = f(x2,w)
    while f1*f2 > 0.0:
        #if x1 >= y:
            #return None,None
        x1 = x2 ; f1 = f2
        x2 = x1+dx; f2 =f(x2,w)
    else:
        return (x1+x2)/2

k=np.linspace(0.01,3.5,100)
for i in range (len(k)):
    s.append(root_stepper(saus,c_Ae,c_A,0.001,k[i]))
for i in range(len(k)):
    kin.append(root_stepper(kink,c_Ae,c_A,0.001,k[i]))
print(kin)
print(s)
plt.plot(k,kin,color='red')
plt.plot(k,s,color='blue')
plt.xlabel('kza')
plt.ylabel('c')
plt.ylim((0,11))
plt.yticks([c_Ae,c_A,],
           ['c_Ae','c_A'])
plt.show()
    
		\end{lstlisting}
\end{spacing}
\begin{figure}
    \flushleft
    \includegraphics[width=6.5cm]{revised fig5.5.png}
    \label{fig:my_label}
\end{figure}
\section{Dispersion Relations under Compressible Effects}
\subsection{Photospheric Condition}
\begin{spacing}
	\begin{lstlisting}[language={Python}]
import matplotlib.pyplot as plt
import numpy as np

x=[]
sol1=[]
sol2=[]
c_s = 2.
c_A = 2.*c_s
c_Ae = c_s/2.
c_se = 3.*c_s/2.
c_t = ((c_s**2*c_A**2)/(c_s**2+c_A**2))**(1/2)
rho_e = 1.0
rho_0 = 1.6*rho_e

def me(c):
    val=(c_se**2-c**2)*(c_Ae**2-c**2)/(c_se**2*c_Ae**2-c**2*(c_se**2+c_Ae**2))
    return np.sqrt(val)

def m0(c):
    val=(c_s**2-c**2)*(c_A**2-c**2)/(c_s**2*c_A**2-c**2*(c_s**2+c_A**2))
    return np.sqrt(abs(val))

def n0(c):
#using n0**2 = -m0**2 
    val=(c**2-c_s**2)*(c**2-c_A**2)/(-c_s**2*c_A**2+c**2*(c_s**2+c_A**2))
    return np.sqrt(abs(val))


def saus_surface(c,kza):
    tmp1 = rho_e*c_Ae**2/me(c) + rho_0*c_A**2/(m0(c)*np.tanh(m0(c)*kza))
    tmp2 = rho_e/me(c)+rho_0/(m0(c)*np.tanh(m0(c)*kza))
    return c**2 - tmp1/tmp2


def kink_surface(c,kza):
    tmp1 = rho_e*c_Ae**2/me(c) + rho_0*c_A**2*np.tanh(m0(c)*kza)/m0(c)
    tmp2 = rho_e/me(c)+rho_0*np.tanh(m0(c)*kza)/m0(c)
    return c**2 - tmp1/tmp2

def saus_body(c,kza):
    tmp1 = rho_e*c_Ae**2*n0(c) + rho_0*c_A**2*me(c)*np.tan(n0(c)*kza)
    tmp2 = rho_e*n0(c) + rho_0*me(c)*np.tan(n0(c)*kza)
    return c**2 - tmp1/tmp2

def kink_body(c,kza):
    tmp1 = rho_e*c_Ae**2/me(c) - rho_0*c_A**2*np.tan(n0(c)*kza)/n0(c)
    tmp2 = rho_e/me(c) - rho_0*np.tan(n0(c)*kza)/n0(c)
    return c**2 - tmp1/tmp2



def root_stepper(f,x,y,dx,k):
    x1 = x ; f1 = f(x,k)
    x2 = x +dx; f2 = f(x2,k)
    while f1*f2 > 0.0:
        #if x1 >= y:
            #return None,None
        x1 = x2 ; f1 = f2
        x2 = x1+dx; f2 =f(x2,k)
    else:
        return (x1+x2)/2
    

s_surface_fast=[]
s_surface_slow=[]
k_surface_fast=[]
k_surface_slow=[]
s_body1=[]
s_body2=[]
k_body=[]
k_body1=[]
k = np.linspace(0.01,3.5,100)

##checking graphs for roots
#cc = np.linspace(c_t,1.85,1000)

#plt.plot(cc,kink_body(cc,0.5))
#plt.plot(cc,kink_body(cc,0.8))
                      
#plt.ylim([-40,40])
#plt.show()
#raise
k1 = np.linspace(0.9,3.5,20)
for i in range (len(k1)):
    k_body.append(root_stepper(kink_body,c_s-0.05,c_t,-0.001,k1[i]))
#print(k_body)
k2 = np.linspace(0.05,0.9,4)
for i in range (len(k2)):
    k_body1.append(root_stepper(kink_body,1.81,c_t,-0.001,k2[i]))
#raise
#print(k_body1)
s1 = np.linspace(0.9,3.5,20)
for i in range (len(s1)):
    s_body1.append(root_stepper(saus_body,c_s-0.01,c_t,-0.001,s1[i]))
s2 = np.linspace(0.05,0.9,5)
for i in range (len(s2)):
    s_body2.append(root_stepper(saus_body,1.85,c_t,-0.001,s2[i]))
    
    
for i in range(len(k)):
    s_surface_fast.append(root_stepper(saus_surface,c_s+0.01,c_se,0.001,k[i]))

for i in range(len(k)):
    s_surface_slow.append(root_stepper(saus_surface,c_Ae+0.1,3,0.001,k[i]))

for i in range(len(k)):
    k_surface_fast.append(root_stepper(kink_surface,c_s-0.19,c_se,0.001,k[i]))

for i in range(len(k)):
    k_surface_slow.append(root_stepper(kink_surface,c_Ae+0.01,c_s,0.001,k[i]))
#for i in range (len(k)):
    #s_body.append(root_stepper(saus_body,c_t,c_s,0.001,k[i]))
#for i in range (len(k)):
    #k_body.append(root_stepper(kink_body,c_t-0.01,c_s,0.001,k[i]))
#set k=0 value
kk=np.zeros(len(k)+1)
kk[0]=0.0
for i in range(1,len(k)+1):
    kk[i] = k[i-1]
    

ss_slow=np.zeros(len(k)+1)
ss_slow[0]=c_t  #from asymptotics
for i in range(1,len(k)+1):
    ss_slow[i] = s_surface_slow[i-1]
   
ss_fast=np.zeros(len(k)+1)
ss_fast[0]=c_se  #from asymptotics
for i in range(1,len(k)+1):
    ss_fast[i] = s_surface_fast[i-1]
    
k_surface_fast = k_surface_fast[6:100]
ks_fast=np.zeros(len(k_surface_fast)+1)
ks_fast[0]=c_t  #from asymptotics
for i in range(1,len(k_surface_fast)+1):
    ks_fast[i] = k_surface_fast[i-1]
kk1 = k[6:100]
ksf = np.zeros(len(k_surface_fast)+1)
for i in range(1,len(k_surface_fast)+1):
    ksf[i] = kk1[i-1]


ks_slow=np.zeros(len(k)+1)
ks_slow[0]=c_Ae  #from asymptotics
for i in range(1,len(k)+1):
    ks_slow[i] = k_surface_slow[i-1]


plt.plot([0,3.5],[c_s,c_s],color='green',linestyle='dashed')
plt.plot([0,3.5],[c_se,c_se],color='green',linestyle='dashed')
plt.plot([0,3.5],[c_Ae,c_Ae],color='green',linestyle='dashed')


plt.ylim([c_Ae-0.2,c_se+0.2])
plt.xlabel('kz*a')
plt.ylabel('c')
plt.yticks([c_t,c_s,c_se,c_Ae],
           ['c_t','c_s','c_se','c_Ae'])


plt.plot(ksf,ks_fast,color='black',linestyle='dashed')
plt.plot(kk,ks_slow,color='black',linestyle='dashed')
plt.plot(kk,ss_fast,color='black')
plt.plot(kk,ss_slow,color='black')
plt.plot(k1,k_body,color='red')
plt.plot(k2,k_body1,color='red')
plt.plot(s1,s_body1,color='blue')
plt.plot(s2,s_body2,color='blue')
plt.show()
#print(s_body)
#print(k_body)

		\end{lstlisting}
\end{spacing}
\begin{figure}
    \flushleft
    \includegraphics[width=6.5cm]{fig5.62.png}
    \label{fig:my_label}
\end{figure}
\subsection{Coronal Condition}
\begin{spacing}
	\begin{lstlisting}[language={Python}]
import matplotlib.pyplot as plt
import numpy as np

rho_0=1
c_se=1/2
c_A=2
c_s=1
c_Ae=5
gamma = 2/11
c_t = ((c_s**2*c_A**2)/(c_s**2+c_A**2))**(1/2)
rho_e = rho_0* ((c_s**2+gamma*(1/2)*c_A**2)/(c_se**2+gamma*(1/2)*c_Ae**2))
#print(rho_e)
#print(c_t)

def m0(c):
    val=(c_s**2-c**2)*(c_A**2-c**2)/(c_s**2*c_A**2-c**2*(c_s**2+c_A**2))
    return np.sqrt(abs(val))

def n0(c):
    val=(c**2-c_s**2)*(c**2-c_A**2)/(-c_s**2*c_A**2+c**2*(c_s**2+c_A**2))
    return np.sqrt(abs(val))
def me(c):
    val= (1-(c**2/c_se**2))
    return np.sqrt(abs(val))


def saus_surf(c,kza):
    a=(c_A/c)**2
    b = 1+ (rho_e/rho_0)*(m0(c)/me(c))*np.tanh(m0(c)*kza)
    return a-b

def kink_surf(c,kza):
    a=me(c)*np.tanh(m0(c)*kza)*(c_A/c)**2
    b= np.tanh(m0(c)*kza)*me(c)+(rho_e/rho_0)*m0(c)
    return a-b

def saus_body(c,kza):
     a=(c_A/c)**2
     b=1-(rho_e/rho_0)*(n0(c)/me(c))*np.tan(n0(c)*kza)
     return a-b

def kink_body(c,kza):
    a=me(c)*np.tan(n0(c)*kza)*(c_A/c)**2
    b=me(c)*np.tan(n0(c)*kza)-(rho_e/rho_0)*n0(c)
    return a-b

def root_stepper(f,x,y,dx,k):
    x1 = x ; f1 = f(x,k)
    x2 = x +dx; f2 = f(x2,k)
    while f1*f2 > 0.0:
        #if x1 >= y:
            #return None,None
        x1 = x2 ; f1 = f2
        x2 = x1+dx; f2 =f(x2,k)
    else:
        return (x1+x2)/2
    
    
k_body_slow=[]
k_body_slow1=[]
s_body_slow=[]
s_body_slow1=[]

k_body_fast=[]
k_body_fast1=[]
k_body_fast2=[]
k_body_fast3=[]
k_body_fast4=[]

s_body_fast=[]
s_body_fast2=[]


k =np.linspace(0.3,3.5,100)

kbs1=np.linspace(0.01,0.3,100)
sbs1=np.linspace(0.01,0.3,100)

kbf=np.linspace(0.01,0.1,100)
kbf1=np.linspace(0.08,3.5,100)
kbf2=np.linspace(0.1,3.5,100)
kbf3=np.linspace(2.3,3.3,100)
kbf4=np.linspace(3.3,3.5,100)

s=np.linspace(0.85,3.5,100)
sbf1=np.linspace(2.2,3.5,100)

for i in range(len(k)):
    k_body_slow.append(root_stepper(kink_body,c_s-0.01,c_t,-0.001,k[i]))
for i in range(len(s)):
    s_body_slow.append(root_stepper(saus_body,c_s-0.01,c_t,-0.001,s[i]))
for i in range (len(kbs1)):
    k_body_slow1.append(root_stepper(kink_body,c_s-0.1,c_t,-0.001,kbs1[i]))
for i in range(len(sbs1)):
    s_body_slow1.append(root_stepper(saus_body,c_s-0.01,c_t,-0.001,sbs1[i]))
#print(k_body_slow1) 

for i in range (len(kbf)):
    k_body_fast.append(root_stepper(kink_body,c_Ae-0.01,4,-0.001,kbf[i]))
for i in range(len(kbf1)):
    k_body_fast1.append(root_stepper(kink_body,c_Ae-0.01,c_A,-0.001,kbf1[i]))
for i in range(len(k)):
    k_body_fast2.append(root_stepper(kink_body,2.2,c_A,-0.001,kbf1[i]))
for i in range (len(kbf3)):
    k_body_fast3.append(root_stepper(kink_body,3.8,c_A,-0.001,kbf3[i]))
for i in range (len(kbf4)):
    k_body_fast4.append(root_stepper(kink_body,3.1,c_A,-0.001,kbf4[i]))


for i in range(len(s)):
    s_body_fast.append(root_stepper(saus_body,c_Ae-0.01,c_A,-0.001,s[i]))
for i in range (len(sbf1)):
    s_body_fast2.append(root_stepper(saus_body,2.8,c_A,-0.001,sbf1[i]))
#print(k_body_fast3)
#print(k_body_fast)

plt.ylim([c_t-0.2,c_Ae+0.2])
plt.xlabel('kz*a')
plt.ylabel('c')
plt.yticks([c_t,c_s,c_Ae,c_A],
           ['c_t','c_s','c_Ae','c_A'])

plt.plot([0,3.5],[c_s,c_s],color='green',linestyle='dashed')
plt.plot([0,3.5],[c_se,c_se],color='green',linestyle='dashed')
plt.plot([0,3.5],[c_Ae,c_Ae],color='green',linestyle='dashed')
plt.plot([0,3.5],[c_A,c_A],color='green',linestyle='dashed')


plt.plot(k,k_body_slow,color='black',linestyle='dashed')
plt.plot(k,s_body_slow,color='black')
plt.plot(kbs1,k_body_slow1,color='black', linestyle='dashed')
plt.plot(sbs1,s_body_slow1,color='black')



plt.plot(kbf,k_body_fast,color='black',linestyle='dashed')
plt.plot(kbf1,k_body_fast1,color='black',linestyle='dashed')
plt.plot(kbf2,k_body_fast2,color='black',linestyle='dashed')
plt.plot(kbf3,k_body_fast3,color='black',linestyle='dashed')
plt.plot(kbf4,k_body_fast4,color='black',linestyle='dashed')


plt.plot(s,s_body_fast,color='black')
plt.plot(sbf1,s_body_fast2,color='black')


plt.show()

		\end{lstlisting}
\end{spacing}
\begin{figure}
    \flushleft
    \includegraphics[width=6.5cm]{coronal_cond.png}
    \label{fig:my_label}
\end{figure}
\subsection{Isolated Slab}
\begin{spacing} 
	\begin{lstlisting}[language={Python}]
import matplotlib.pyplot as plt
import numpy as np

x=[]
sol1=[]
sol2=[]
c_s = 2.
c_A = 2.*c_s
c_Ae = 0
c_se = 3.*c_s/2.
c_t = ((c_s**2*c_A**2)/(c_s**2+c_A**2))**(1/2)
rho_e = 1.0
rho_0 = 1.6*rho_e

def me(c):
    val=(c_se**2-c**2)*(c_Ae**2-c**2)/(c_se**2*c_Ae**2-c**2*(c_se**2+c_Ae**2))
    return np.sqrt(val)

def m0(c):
    val=(c_s**2-c**2)*(c_A**2-c**2)/(c_s**2*c_A**2-c**2*(c_s**2+c_A**2))
    return np.sqrt(abs(val))

def n0(c):
#using n0**2 = -m0**2 
    val=(c**2-c_s**2)*(c**2-c_A**2)/(-c_s**2*c_A**2+c**2*(c_s**2+c_A**2))
    return np.sqrt(abs(val))


def saus_surface(c,kza):
    tmp1 = rho_e*c_Ae**2/me(c) + rho_0*c_A**2/(m0(c)*np.tanh(m0(c)*kza))
    tmp2 = rho_e/me(c)+rho_0/(m0(c)*np.tanh(m0(c)*kza))
    return c**2 - tmp1/tmp2


def kink_surface(c,kza):
    tmp1 = rho_e*c_Ae**2/me(c) + rho_0*c_A**2*np.tanh(m0(c)*kza)/m0(c)
    tmp2 = rho_e/me(c)+rho_0*np.tanh(m0(c)*kza)/m0(c)
    return c**2 - tmp1/tmp2

def saus_body(c,kza):
    tmp1 = rho_e*c_Ae**2*n0(c) + rho_0*c_A**2*me(c)*np.tan(n0(c)*kza)
    tmp2 = rho_e*n0(c) + rho_0*me(c)*np.tan(n0(c)*kza)
    return c**2 - tmp1/tmp2

def kink_body(c,kza):
    tmp1 = rho_e*c_Ae**2/me(c) - rho_0*c_A**2*np.tan(n0(c)*kza)/n0(c)
    tmp2 = rho_e/me(c) - rho_0*np.tan(n0(c)*kza)/n0(c)
    return c**2 - tmp1/tmp2



def root_stepper(f,x,y,dx,k):
    x1 = x ; f1 = f(x,k)
    x2 = x +dx; f2 = f(x2,k)
    while f1*f2 > 0.0:
        #if x1 >= y:
            #return None,None
        x1 = x2 ; f1 = f2
        x2 = x1+dx; f2 =f(x2,k)
    else:
        return (x1+x2)/2
    

s_surface_fast=[]
s_surface_slow=[]
k_surface_fast=[]
k_surface_slow=[]
s_body1=[]
s_body2=[]
k_body=[]
k_body1=[]
k = np.linspace(0.01,3.5,100)

#cc = np.linspace(c_t,1.85,1000)

#plt.plot(cc,kink_body(cc,0.5))
#plt.plot(cc,kink_body(cc,0.8))
                      
#plt.ylim([-40,40])
#plt.show()
#raise
k1 = np.linspace(0.9,3.5,20)
for i in range (len(k1)):
    k_body.append(root_stepper(kink_body,c_s-0.05,c_t,-0.001,k1[i]))
#print(k_body)
k2 = np.linspace(0.05,0.9,4)
for i in range (len(k2)):
    k_body1.append(root_stepper(kink_body,1.81,c_t,-0.001,k2[i]))
#raise
#print(k_body1)
s1 = np.linspace(0.9,3.5,20)
for i in range (len(s1)):
    s_body1.append(root_stepper(saus_body,c_s-0.01,c_t,-0.001,s1[i]))
s2 = np.linspace(0.05,0.9,5)
for i in range (len(s2)):
    s_body2.append(root_stepper(saus_body,1.85,c_t,-0.001,s2[i]))
    
    
for i in range(len(k)):
    s_surface_fast.append(root_stepper(saus_surface,c_s+0.01,c_se,0.001,k[i]))

for i in range(len(k)):
    s_surface_slow.append(root_stepper(saus_surface,c_Ae+0.1,3,0.001,k[i]))

for i in range(len(k)):
    k_surface_fast.append(root_stepper(kink_surface,c_s-0.19,c_se,0.001,k[i]))

for i in range(len(k)):
    k_surface_slow.append(root_stepper(kink_surface,c_Ae+0.01,c_s,0.001,k[i]))
#for i in range (len(k)):
    #s_body.append(root_stepper(saus_body,c_t,c_s,0.001,k[i]))
#for i in range (len(k)):
    #k_body.append(root_stepper(kink_body,c_t-0.01,c_s,0.001,k[i]))

kk=np.zeros(len(k)+1)
kk[0]=0.0
for i in range(1,len(k)+1):
    kk[i] = k[i-1]
    

ss_slow=np.zeros(len(k)+1)
ss_slow[0]=c_t  #from asymptotics
for i in range(1,len(k)+1):
    ss_slow[i] = s_surface_slow[i-1]
   
ss_fast=np.zeros(len(k)+1)
ss_fast[0]=c_se  #from asymptotics
for i in range(1,len(k)+1):
    ss_fast[i] = s_surface_fast[i-1]
    
k_surface_fast = k_surface_fast[6:100]
ks_fast=np.zeros(len(k_surface_fast)+1)
ks_fast[0]=c_t  
for i in range(1,len(k_surface_fast)+1):
    ks_fast[i] = k_surface_fast[i-1]
#make k just for this
kk1 = k[6:100]
ksf = np.zeros(len(k_surface_fast)+1)
for i in range(1,len(k_surface_fast)+1):
    ksf[i] = kk1[i-1]


ks_slow=np.zeros(len(k)+1)
ks_slow[0]=c_Ae  #from asymptotics
for i in range(1,len(k)+1):
    ks_slow[i] = k_surface_slow[i-1]

plt.plot([0,3.5],[c_s,c_s],color='green',linestyle='dashed')
plt.plot([0,3.5],[c_se,c_se],color='green',linestyle='dashed')
plt.plot([0,3.5],[c_Ae,c_Ae],color='green',linestyle='dashed')


plt.ylim([c_Ae-0.2,c_se+0.2])
plt.xlabel('kz*a')
plt.ylabel('c')
plt.yticks([c_t,c_s,c_se,c_Ae],
           ['c_t','c_s','c_se','c_Ae'])


plt.plot(ksf,ks_fast,color='black',linestyle='dashed')
plt.plot(kk,ks_slow,color='black',linestyle='dashed')
plt.plot(kk,ss_fast,color='black')
plt.plot(kk,ss_slow,color='black')
plt.plot(k1,k_body,color='red')
plt.plot(k2,k_body1,color='red')
plt.plot(s1,s_body1,color='blue')
plt.plot(s2,s_body2,color='blue')
plt.show()
#print(s_body)
#print(k_body)
