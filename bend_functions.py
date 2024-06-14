import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
import numpy as np
from scipy.special import fresnel
from scipy import special

def euler_bend(R, theta, num_points):
    L = 2 * abs(R) * theta # total length of the Euler bend
    t = np.linspace(0, L, num_points//2)
    f = np.sqrt(np.pi * abs(R) * L) # Fresnel integral function are defined as function of (pi*t^2/2)
    y, x = fresnel(t / f) 
    print('Total lenght of the Euler:',L, 'um')
    return f*x, f*y


def euler_bend_NA(R, theta, N, num_points):
    x_data=np.array([])
    y_data=np.array([])
    L = (N+1) * R * theta
    s = np.linspace(0, L, num_points//2)
    k=theta*(N+1)/L**(N+1)
    #f = np.sqrt(np.pi * R * L**(N)*(N+1)/2) + 1e-18 # for numerical stability
    for j in range(0,len(s)):
        x=integrate.quad(lambda s: np.cos(k*s**(N+1)/(N+1)), 0, s[j])[0]
        x_data=np.append(x_data,x)
        y=integrate.quad(lambda s: np.sin(k*s**(N+1)/(N+1)), 0, s[j])[0]
        y_data=np.append(y_data,y)
        
    print('Total lenght of the NA-Euler with n=',N,',:',L, 'um')
    return x_data, y_data




def euler_bend_NA_90(Reff, n, num_points):
    k=np.linspace(1e-8*10**(-5*(n+0.1)),0.02*Reff,100000)
    x_data=np.array([])
    y_data=np.array([])

    R_min=Reff-0.01*np.abs(n+0.1)
    R_max=Reff+0.01*np.abs(n+0.1)
    min_result = 0
    k1 = 0

    for i in range(0,len(k)):
        x=integrate.quad(lambda l: np.cos(k[i]*l**(n+1)/(n+1)), 0, ((np.pi*(n+1)/(4*k[i])))**(1/(n+1)))[0]
        y=integrate.quad(lambda l: np.sin(k[i]*l**(n+1)/(n+1)), 0, ((np.pi*(n+1)/(4*k[i])))**(1/(n+1)))[0]
        result=x+y
        
        if result < R_max and result > R_min:
            k1=k[i]
        elif result < R_min and result > min_result:
            kmax=k[i]
            min_result=result
        elif result > R_max:
            kmin=k[i]

    if k1 == 0:
        k=np.linspace(kmin,kmax,10000)
        for i in range(0,len(k)):
            R_min=Reff-0.01*np.abs(n+0.1)
            R_max=Reff+0.01*np.abs(n+0.1)
            x=integrate.quad(lambda l: np.cos(k[i]*l**(n+1)/(n+1)), 0, ((np.pi*(n+1)/(4*k[i])))**(1/(n+1)))[0]
            y=integrate.quad(lambda l: np.sin(k[i]*l**(n+1)/(n+1)), 0, ((np.pi*(n+1)/(4*k[i])))**(1/(n+1)))[0]
            result=x+y

            if result < R_max and result > R_min:
                k1=k[i]
            elif result < R_min:
                kmax=k[i]
            elif result > R_max:
                kmin=k[i]

    
    s_f=((np.pi*(n+1)/(4*k1)))**(1/(n+1))
    print('Total lenght of the 90-Euler with n=',n,',:',2*s_f, 'um')
    s_p=np.linspace(0, s_f, num_points//2)

    for j in range(0,len(s_p)):
        x=integrate.quad(lambda s: np.cos(k1*s**(n+1)/(n+1)), 0, s_p[j])[0]
        x_data=np.append(x_data,x)
        y=integrate.quad(lambda s: np.sin(k1*s**(n+1)/(n+1)), 0, s_p[j])[0]
        y_data=np.append(y_data,y)
    
    x2_data, y2_data = np.dot(np.array([[0, 1], [-1, 0]]), np.stack([x_data, y_data], 0))

    x2_data, y2_data = -x2_data[::-1], y2_data[::-1]

    x2_data, y2_data = x2_data-x2_data[0]+x_data[-1], y2_data-y2_data[0]+y_data[-1]

    x_data1 = np.concatenate([x_data, x2_data], 0)
    y_data1 = np.concatenate([y_data, y2_data], 0)
    
    return x_data1, y_data1, s_f



def euler_bend_NA_90A(Reff, n, num_points): # old version. It's faster and it works, but it's not completely correct from the mathematical point of view (ref. [2])
    a=0
    k=1
    A=np.arange(0.5,Reff*2,0.01)
    s_p=np.linspace(0, ((np.pi*(n+1)/(4*k)))**(1/(n+1)), num_points//2)
    x_data=np.array([])
    y_data=np.array([])

    for i in range(0,len(A)):
        b=(np.pi*(n+1)/4)**(1/(n+1))
        R_min=Reff-0.01
        R_max=Reff+0.01
        x=integrate.quad(lambda l: np.cos(l**(n+1)/(n+1)), a, b)[0]
        y=integrate.quad(lambda l: np.sin(l**(n+1)/(n+1)), a, b)[0]
        result=A[i]*x+A[i]*y
        #print(A[i])
        #print(x,y)
        #print(result)
        if result < R_max and result > R_min:
            #print("Constant A=", A[i])
            A1=A[i]
    for j in range(0,len(s_p)):
        x=integrate.quad(lambda s: np.cos(k*s**(n+1)/(n+1)), a, s_p[j])[0]
        x_data=np.append(x_data,x*A1)
        y=integrate.quad(lambda s: np.sin(k*s**(n+1)/(n+1)), a, s_p[j])[0]
        y_data=np.append(y_data,y*A1)
    
    # first, rotate by the final angle
    x2_data, y2_data = np.dot(np.array([[0, 1], [-1, 0]]), np.stack([x_data, y_data], 0))
    # then, flip along the x-axis (and reverse direction of the curve):
    x2_data, y2_data = -x2_data[::-1], y2_data[::-1]
    # then translate from (x2[0], y2[0]) to (x1[-1], y1[-1])
    x2_data, y2_data = x2_data-x2_data[0]+x_data[-1], y2_data-y2_data[0]+y_data[-1]

    x_data = np.concatenate([x_data, x2_data], 0)
    y_data = np.concatenate([y_data, y2_data], 0)
    

    return x_data,y_data



#
#=============Cubic Bezier curve====================
#


class Point(object):
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y


class CubicBezier(object):
    def __init__(self, p0x= 0, p0y= 0, p1x= 0, p1y= 0, p2x= 0, p2y= 0, p3x= 0, p3y= 0):
        self.p0 = Point(p0x, p0y)
        self.p1 = Point(p1x, p1y)
        self.p2 = Point(p2x, p2y)
        self.p3 = Point(p3x, p3y)
        self.obstacles = []

    def max_k(self, granuality=100):
        'Calculate maximal curvature of the cubic Bezier curve.'
        k = 0
        for t in range(0, granuality):
            t = t / granuality
            x_d = 3 * ((1 - t) ** 2) * (self.p1.x - self.p0.x) + 6 * (1 - t) * t * (self.p2.x - self.p1.x) + 3 * (t ** 2) * (
                        self.p3.x - self.p2.x)
            y_d = 3 * ((1 - t) ** 2) * (self.p1.y - self.p0.y) + 6 * (1 - t) * t * (self.p2.y - self.p1.y) + 3 * (t ** 2) * (
                        self.p3.y - self.p2.y)
            x_dd = 6 * (1 - t) * (self.p2.x - 2 * self.p1.x + self.p0.x) + 6 * t * (self.p3.x - 2 * self.p2.x + self.p1.x)
            y_dd = 6 * (1 - t) * (self.p2.y - 2 * self.p1.y + self.p0.y) + 6 * t * (self.p3.y - 2 * self.p2.y + self.p1.y)
            k = max(k,abs(x_d*y_dd - y_d*x_dd)/math.pow(x_d**2 + y_d**2, 3/2))
        return k

    def calc_curve(self, granuality=100):
        'Calculate the cubic Bezier curve with the given granuality.'
        B_x = []
        B_y = []
        for t in range(0, granuality):
            t = t / granuality
            x = ((1 - t) ** 3) * self.p0.x + 3 * ((1 - t) ** 2) * t * self.p1.x + 3 * (1 - t) * (t ** 2) * self.p2.x\
                + (t ** 3) * self.p3.x
            y = ((1 - t) ** 3) * self.p0.y + 3 * ((1 - t) ** 2) * t * self.p1.y + 3 * (1 - t) * (t ** 2) * self.p2.y\
                + (t ** 3) * self.p3.y
            B_x.append(x)
            B_y.append(y)
        return [B_x, B_y]



############################ 180 ###########################
def euler_bend_Hybrid_180(Reff, p, num_points): 
    a=0
    k=(2*np.sqrt(2)*(integrate.quad(lambda t: np.sin(t**2), 0, np.sqrt(np.pi*p/2))[0])+2/np.sqrt(np.pi*p)*np.sin(np.pi*(1-p)/2))**2/(2*Reff)**2
    #Rf=(k*np.pi*p)**(-1)
    #s_p=np.linspace(0, (Rf*k)**(-1/2), num_points//2)
    
    s_p=np.linspace(0, np.sqrt(np.pi*p/k), num_points//2)
    x_data=np.array([])
    y_data=np.array([])

    for j in range(0,len(s_p)):
        x=integrate.quad(lambda l: np.cos(k*l**2/(2)), a, s_p[j])[0]
        x_data=np.append(x_data,x)
        y=integrate.quad(lambda l: np.sin(k*l**2/(2)), a, s_p[j])[0]
        y_data=np.append(y_data,y)
    
    x0=x_data[-1]-np.sin(np.pi*p/2)/np.sqrt(k*np.pi*p)
    # first, rotate by the final angle
    x2_data, y2_data = np.dot(np.array([[-1, 0], [0, -1]]), np.stack([x_data, y_data], 0))
    # then, flip along the x-axis (and reverse direction of the curve):
    x2_data, y2_data = -x2_data[::-1], y2_data[::-1]
    # then translate from (x2[0], y2[0]) to (x1[-1], y1[-1])
    x2_data, y2_data = x2_data, y2_data + 2*Reff

    theta = np.linspace( -(-p+1)*np.pi/2, (-p+1)*np.pi/2 , 150 )
    radius = 1/np.sqrt(k*np.pi*p)
    a = (radius * np.cos( theta ))+x0
    b = (radius * np.sin( theta ))+Reff

    x_data = np.concatenate([x_data, a,x2_data], 0)
    y_data = np.concatenate([y_data, b,y2_data], 0)

    s_f= 2*np.sqrt(np.pi*p/k)+radius*(-p+1)*np.pi
    print('Total lenght of the Hybrid-Euler with p=',p,',:',s_f, 'um')
    return x_data,y_data

def euler_bend_NA_180(Reff, n, num_points):
    k=np.linspace(1e-8*1**(-5*(n+0.1)),0.02*Reff,100)
    x_data=np.array([])
    y_data=np.array([])

    R_min=Reff-0.01*np.abs(n+0.1)
    R_max=Reff+0.01*np.abs(n+0.1)
    min_result = 0
    k1 = 0
    kmin = 1e-18

    while k1 == 0:
      for i in range(0,len(k)):
        x=integrate.quad(lambda l: np.cos(k[i]*l**(n+1)/(n+1)), 0, ((np.pi*(n+1)/(4*k[i])))**(1/(n+1)))[0]
        y=integrate.quad(lambda l: np.sin(k[i]*l**(n+1)/(n+1)), 0, ((np.pi*(n+1)/(4*k[i])))**(1/(n+1)))[0]
        result=x+y
        
        if result < R_max and result > R_min:
            k1=k[i]
        elif result < R_min and result > min_result:
            kmax=k[i]
            min_result=result
        elif result > R_max:
            kmin=k[i]
      k=np.linspace(kmin,kmax,10000)
    
    s_f=((np.pi*(n+1)/(4*k1)))**(1/(n+1))
    s_p=np.linspace(0, s_f, num_points//2)

    for j in range(0,len(s_p)):
        x=integrate.quad(lambda s: np.cos(k1*s**(n+1)/(n+1)), 0, s_p[j])[0]
        x_data=np.append(x_data,x)
        y=integrate.quad(lambda s: np.sin(k1*s**(n+1)/(n+1)), 0, s_p[j])[0]
        y_data=np.append(y_data,y)
    
    # first, rotate by the final angle
    x2_data, y2_data = np.dot(np.array([[0, 1], [-1, 0]]), np.stack([x_data, y_data], 0))
    # then, flip along the x-axis (and reverse direction of the curve):
    x2_data, y2_data = -x2_data[::-1], y2_data[::-1]
    # then translate from (x2[0], y2[0]) to (x1[-1], y1[-1])
    x2_data, y2_data = x2_data-x2_data[0]+x_data[-1], y2_data-y2_data[0]+y_data[-1]

    x_data1 = np.concatenate([x_data, x2_data], 0)
    y_data1 = np.concatenate([y_data, y2_data], 0)
    

      # first, rotate by the final angle
    x2_data, y2_data = np.dot(np.array([[1, 0], [0, -1]]), np.stack([x_data1, y_data1], 0))
    # then, flip along the x-axis (and reverse direction of the curve):
    x2_data, y2_data = x2_data[::-1], y2_data[::-1]
    # then translate from (x2[0], y2[0]) to (x1[-1], y1[-1])
    x2_data, y2_data = x2_data-x2_data[0]+x_data1[-1], y2_data-y2_data[0]+y_data1[-1]

    x_data = np.concatenate([x_data1, x2_data], 0)
    y_data = np.concatenate([y_data1, y2_data], 0)
    
    print('Total lenght of the 180-Euler with n=',n,',:',4*s_f, 'um')
    return x_data, y_data, 4*s_f


###################### s-bend ############################

def tanh(T,d,L,npts):
    x = np.linspace( -L/2, L/2 , npts )
    y=d/2*np.tanh(x/T)
    #s_f=integrate.quad(lambda x:np.sqrt(1+d**2/(2*T)**2*1/np.cosh(x/T)**4), 0, L)[0]
    #R_av=s_f/np.arctan(d/(2*T))
    #print("The average radius of curvature for the tanh function is {0:.3f} mm".format( R_av*1e-3))
    return x+L/2,y+d/2

def erf(E,d,L,npts):
    x = np.linspace(-L/2, L/2,npts)
    y=d/2*special.erf(x/E)
    #s_f=integrate.quad(lambda x:np.sqrt(1+(d/(np.sqrt(np.pi)*E)*np.exp(-(x/E)**2))**2), 0, L)[0]
    #R_av=s_f/np.arctan(d/(np.sqrt(np.pi)*E))
    #print("The average radius of curvature for the erf function is {0:.3f} mm".format( R_av*1e-3))
    return x+L/2,y+d/2

def arctan(C,d,L,npts):
    x = np.linspace(-L/2, L/2,npts)
    y=2/(np.pi)*np.arctan(np.pi*(x*C)/2)
    #s_f=integrate.quad(lambda x:np.sqrt(1+(d/(C*(1+(np.pi*x/(2*C))**2))**2)), 0, L)[0]
    #R_av=s_f/np.arctan(d/C)
    #print("The average radius of curvature for the arctan function is {0:.3f} mm".format( R_av*1e-3))
    return x+L/2,y*d/(y[npts-1]-y[0])+d/2

def rad(C,d,L,npts):
    x = np.linspace(-L/2, L/2,npts)
    y=x*C/(np.sqrt(1+(x*C)**2))
    #s_f=integrate.quad(lambda x:np.sqrt(1+(d/(C*((x/C)**2+1)**(3/2)))**2), 0, L)[0]
    #R_av=s_f/np.arctan(d/C)
    #print("The average radius of curvature for the radical function is {0:.3f} mm".format( R_av*1e-3))
    return x+L/2,y*d/(y[npts-1]-y[0])+d/2

def abv(C,d,L,npts):
    x = np.linspace(-L/2, L/2,npts)
    y=x*C*1/(1+abs(x*C))
    s_f=integrate.quad(lambda x:np.sqrt(1+(d/C*(-x**2/(C**2*np.abs(x/C))+np.abs(x/C)+1)/(np.abs(x/C)+1)**2)**2), 0, L)[0]
    #R_av=s_f/np.arctan(d/C)
    #print("The average radius of curvature for the absolute value function is {0:.3f} mm".format( R_av*1e-3))
    return x+L/2,y*d/(y[npts-1]-y[0])+d/2


def euler_bend_S2(theta, R, d,  n,num_points):
    #R=d/20
    k= 1/(theta**n*R**(n+1)*(n+1)**n)    #k=1/(2*R**2*theta)
    s_f = (theta*(n+1)/k )**(1/(n+1))    #(np.sqrt(2*theta/k))
    s = np.linspace(0, s_f, num_points//2)

    #t = np.sqrt(np.pi/(k)) 
    #y1, x1 = fresnel(s/ t)
    #y1=t*y1
    #x1=t*x1

    x1=np.array([])
    y1=np.array([])

    for j in range(0,len(s)):
        x=integrate.quad(lambda l: np.cos(k*l**(n+1)/(n+1)), 0, s[j])[0]
        x1=np.append(x1,x)
        y=integrate.quad(lambda l: np.sin(k*l**(n+1)/(n+1)), 0, s[j])[0]
        y1=np.append(y1,y)
    

    # first, rotate by the final angle
    
    x2_data, y2_data = np.dot(np.array([[-1, 0], [0, -1]]), np.stack([x1, y1], 0))
    x2_data, y2_data = x2_data[::-1], y2_data[::-1]
    x2_data, y2_data = x2_data+2*x1[-1], y2_data+2*y1[-1]

    x_data = np.concatenate([x1, x2_data], 0)
    y_data = np.concatenate([y1, y2_data], 0)
    
    L=d*x_data[-1]/y_data[-1]
    x_data=x_data*L/x_data[-1]
    y_data=y_data*d/y_data[-1]
    #R_av=s_f**(-n)/(k*np.abs(1-n)) #s_f*(n+1)/(k*s_f**(n+1))   #(2*theta**(n-1)*(n+1)**(n-1)*R**n)    #(n+1)/(4*k*s_f) #1/(2*k*s_f)
    #print("The average radius of curvature for the 2-Euler curve is {0:.0f} mm, with n= {1},bending angle $\Theta=\pi$/{2:.0f} and radius at the central bending point $R_f=${3} mm".format( R_av*1e-3,n,np.pi/theta,R*1e-3))

    return x_data, y_data


def euler_bend_S4(theta, R, d, n, num_points):
    #R=d/20
    k= 1/(theta**n*R**(n+1)*(n+1)**n)    #1/(2*R**2*theta)
    s_f = (theta*(n+1)/k )**(1/(n+1))    #(np.sqrt(2*theta/k))
    s = np.linspace(0, s_f, num_points)


    x1=np.array([])
    y1=np.array([])

    for j in range(0,len(s)):
        x=integrate.quad(lambda l: np.cos(k*l**(n+1)/(n+1)), 0, s[j])[0]
        x1=np.append(x1,x)
        y=integrate.quad(lambda l: np.sin(k*l**(n+1)/(n+1)), 0, s[j])[0]
        y1=np.append(y1,y)
    
    
    #y1, x1 = fresnel(s/ t)
    #y1=t*y1
    #x1=t*x1
    # first, rotate by the final angle
    
    x2, y2 = np.dot(np.array([[np.cos(2*theta), np.sin(2*theta)], [-np.sin(2*theta), np.cos(2*theta)]]), np.stack([x1, y1], 0))
    # then, flip along the x-axis (and reverse direction of the curve):
    x2, y2 = -x2[::-1], y2[::-1]
    # then translate from (x2[0], y2[0]) to (x1[-1], y1[-1])
    x2, y2 = x2-x2[0]+x1[-1], y2-y2[0]+y1[-1]
    x21 = np.concatenate([x1, x2], 0)
    y21 = np.concatenate([y1, y2], 0)

    x2, y2 = np.dot(np.array([[-1, 0], [0, -1]]), np.stack([x21, y21], 0))
    # then, flip along the x-axis (and reverse direction of the curve):
    x2, y2 = x2[::-1], y2[::-1]
    # then translate from (x2[0], y2[0]) to (x1[-1], y1[-1])
    x2, y2 = x2- x2[0]+x21[-1], y2-y2[0]+y21[-1]
    x = np.concatenate([x21, x2], 0)
    y = np.concatenate([y21, y2], 0)
    

    #R_av=s_f**n/(k*np.abs(1-n)) #s_f*(n+1)/(k*s_f**(n+1))   #(2*theta**(n-1)*(n+1)**(n-1)*R**n)    #(n+1)/(4*k*s_f) #1/(2*k*s_f)
    #print("The average radius of curvature for the 4-Euler curve is {0:.0f} mm, with n= {1},bending angle $\Theta=\pi$/{2:.0f} and radius at the central bending point $R_f=${3} mm".format( R_av*1e-3,n,np.pi/theta,R*1e-3))

    L1=d*x[-1]/y[-1]
    x=x*L1/x[-1]
    y=y*d/y[-1]

    return x, y

def circular(r,t1,t2,x0,y0,npts):
    theta = np.linspace( t1, t2 , npts)
    a = (r * np.cos( theta ))+x0
    b = (r * np.sin( theta ))+y0

def circ_S(r,d,npts):
    thetaf=np.arccos(1-d/(2*r))
    theta = np.linspace( 3*np.pi/2,3*np.pi/2+thetaf, npts//2) #if we want to reach higher values of the slope at the center, we should increase the number of points
    a = (r * np.cos( theta ))
    b = (r * np.sin( theta ))+r

    x2, y2 = np.dot(np.array([[-1, 0], [0, -1]]), np.stack([a, b], 0))
    x2, y2 = x2[::-1], y2[::-1]
    x2, y2 = x2+2*a[-1], y2+2*b[-1]

    a2 = np.concatenate([a, x2], 0)
    b2 = np.concatenate([b, y2], 0)
    return a2,b2




def Series_S2Euler(L,theta, R, d, n,npts):
    x1, y1= euler_bend_S2(theta, R, d, n,npts)
    n=L//x1[-1]
    print("Number of bends (2-Euler):", n, "for L=", L*1e-4 ,"cm")
    y=np.array([])
    x=np.array([])

    for i in range(1,int(n)):
        x2=x1+x1[-1]*i
        y2=((-1)**i)*y1+(i*y1[-1]-((-1)**i)*i*y1[-1])/(2*i)
        y=np.append(y,y2)
        x=np.append(x,x2)

    x=np.append(x1,x)
    y=np.append(y1,y)
    return x,y


def Series_S4Euler(L,theta, R, d, n, npts):
    x1, y1= euler_bend_S4(theta, R, d, n, npts)
    n=L//x1[-1]
    print("Number of S-bends (4-Euler):", n, "for L=", L*1e-4 ,"cm")
    y=np.array([])
    x=np.array([])

    for i in range(1,int(n)):
        x2=x1+x1[-1]*i
        y2=((-1)**i)*y1+(i*y1[-1]-((-1)**i)*i*y1[-1])/(2*i)
        y=np.append(y,y2)
        x=np.append(x,x2)

    x=np.append(x1,x)
    y=np.append(y1,y)
    return x,y 



def Series(curve, Ltot):
    x1, y1= curve
    n=Ltot//x1[-1]
    print("Number of S-bends:", n, "for L=", Ltot*1e-4 ,"cm")
    y=np.array([])
    x=np.array([])

    for i in range(1,int(n)):
        x2=x1+x1[-1]*i
        y2=((-1)**i)*y1+(i*y1[-1]-((-1)**i)*i*y1[-1])/(2*i)
        y=np.append(y,y2)
        x=np.append(x,x2)

    x=np.append(x1,x)
    y=np.append(y1,y)
    return x,y 