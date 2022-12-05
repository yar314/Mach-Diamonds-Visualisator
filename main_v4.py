import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math as mt
from matplotlib.colors import LinearSegmentedColormap
import tkinter as tk

def om(M):
    y=1.4
    return mt.sqrt((y+1)/(y-1))*mt.atan(mt.sqrt((M**2-1)*(y-1)/(y+1)))-mt.atan(mt.sqrt(M**2-1))

def Mach_from_omega(omega):
    omega = omega*(180/3.1415926)
    return 0.000014*mt.pow(omega,3)-0.0011*omega**2+0.0631*omega+0.9447


def pr_ratio(M):
    y=1.4
    return (1+M**2*(y-1)/2)**(-y/(y-1))

def Mach(pr_ratio):
    y=1.4
    return mt.sqrt((pr_ratio**((1-y)/y)-1)*2/(y-1))

def I(t,omega):
    return (t+omega)/2+500

def II(t,omega):
    return (t-omega)/2+500

def theta(I,II):
    return I+II-1000

def omega(I,II):
    return I-II

def perp( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def seg_int(a, b) :
    a=np.array(a)
    b=np.array(b)
    da = a[1]-a[0]
    db = b[1]-b[0]
    dp = a[0]-b[0]
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b[0]

def intersect_3(a,b,c):
    x1 = seg_int(c,b)[0]
    y1 = seg_int(c,b)[1]
    x2 = seg_int(a,c)[0]
    y2 = seg_int(a,c)[1]
    x3 = seg_int(a,b)[0]
    y3 = seg_int(a,b)[1]
    return [[x1,y1],[x2,y2],[x3,y3]]

def intersect_4(a,b,c,d):
    x1 = seg_int(a,b)[0]
    y1 = seg_int(a,b)[1]
    x2 = seg_int(b,c)[0]
    y2 = seg_int(b,c)[1]
    x3 = seg_int(c,d)[0]
    y3 = seg_int(c,d)[1]
    x4 = seg_int(d,a)[0]
    y4 = seg_int(d,a)[1]
    return [[x1,y1],[x2,y2],[x3,y3],[x4,y4]]

def start():
    plt.axis('off')  # command for hiding the axis.
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    #Data-----------------------------------------------------------------------------------------------------------------
    M=float(entry_M.get())
    ratio=float(entry_ratio.get())
    n=int(entry_n.get())
    y=1.4
    t_0=0
    omega_0 = om(M)
    pr_ratio_0= pr_ratio(M)
    I_0=I(t_0,omega_0)
    II_0=II(t_0,omega_0)
    pr_ratio_n=ratio*pr_ratio_0
    M_n = Mach(pr_ratio_n)
    omega_n = om(M_n)
    t_n=omega_0-omega_n
    I_n=I_0
    II_n=II(t_n,omega_n)

    if n%2 ==0:
        theta=[t_n]
        for i in range(int((n-1)/2)):
            theta.append(2*(I(theta[i],omega_n)-theta[i]/(n-i))-omega_n-1000)
        theta.append(0)
        for i in range(int((n-1)/2)-1,-1,-1):
            theta.append(-(2*(I(theta[i],omega_n)-theta[i]/(n-i))-omega_n-1000))
        theta.append(-t_n)
    else:
        theta=[t_n]
        for i in range(int((n-1)/2)):
            theta.append(2*(I(theta[i],omega_n)-theta[i]/(n-i))-omega_n-1000)
        for i in range(int((n-1)/2)-1,-1,-1):
            theta.append(-(2*(I(theta[i],omega_n)-theta[i]/(n-i))-omega_n-1000))
        theta.append(-t_n)


    #Slip boundary-------------------------------------------------------------------------------------------------------------
    x=[10]
    y1=[0]
    y2=[80]
    k=120/n
    for i in range(len(theta)):
        s_x=np.array([x[i],mt.cos(theta[i])*k+x[i]])
        x.append(mt.cos(theta[i])*k+x[i])
        s1_y=np.array([y1[i],mt.sin(theta[i])*k+y1[i]])
        y1.append(mt.sin(theta[i])*k+y1[i])
        s2_y=np.array([y2[i],-mt.sin(theta[i])*k+y2[i]])
        y2.append(-mt.sin(theta[i])*k+y2[i])
        plt.plot(s_x,s1_y,'k--',lw=1)
        plt.plot(s_x,s2_y,'k--',lw=1)

    def distribution():
        #Color-----------------------------------------------------------------------------------------------------------------------
        m_grad = []
        for i in range(n):
            m=Mach_from_omega(2*I_0-1000-(t_n/n)*i)
            m_grad.append(m)
        for i in range(len(x)-3):
            m=Mach_from_omega(omega_n)
            m_grad.append(m)
        for i in range(n-1):
            for j in range(n-1-i):
                m=Mach_from_omega(omega(I(theta[i+1],omega_n),(j+2)*t_n/n+1000-I(theta[i],omega_n)))
                m_grad.append(m)
            m=Mach_from_omega(omega(I(theta[i+1],omega_n),t_n+1000-I(t_n,omega_n)))
            m_grad.append(m)
            for j in range(i):
                m=Mach_from_omega(omega(I(theta[i+1],omega_n),theta[1+j]+1000-I(theta[j+1],omega_n)))
                m_grad.append(m)
        m=Mach_from_omega(omega(I(-theta[0],omega_n),II(theta[0],omega_n)))
        m_grad.append(m)

        m_max=mt.ceil(max(m_grad)*100)/100
        m_min=mt.floor(min(m_grad)*100)/100

        r='ff'
        g='ff'
        b='99'

        k0=16+(int(g,16)-16)*m_min/(m_min-m_max)
        k1=-(int(g,16)-16)/(m_min-m_max)
        #Colorbar------------------------------------------------------------------------------------------------------------

        colors = [(int(r,16)/255, 16/255, int(b,16)/255),(int(r,16)/255, int(g,16)/255, int(b,16)/255)]  # RGB

        cmap = LinearSegmentedColormap.from_list('mach_colors', colors, N=100)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(m_min, m_max))
        plt.colorbar(sm,orientation='vertical',label="Mach number",fraction=0.03,location='left',aspect=15,pad=0.01)

        #Filling----------------------------------------------------------------------------------------------------------------------
        dots_low = []
        dots_up = []
        for i in range(len(x)):
            dots_low.append([x[i],y1[i]])
            dots_up.append([x[i],y2[i]]) 

        for i in range(n):
            m=Mach_from_omega(2*I_0-1000-(t_n/n)*i)
            color='#'+r+str(hex(int(k1*m+k0))).split('x')[-1]+b
            if i==0:
                plt.gca().add_patch(plt.Polygon([[-20,80],[-20,0],[10,0],[10,80]],color=color))
                plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[i],dots_up[0]],[dots_up[1],dots_low[0]],[dots_up[0],dots_low[i+1]]),color=color))
                plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[len(x)-1],dots_up[len(x)-1-i]],[dots_up[len(x)-1],dots_low[len(x)-2]],[dots_up[len(x)-2-i],dots_low[len(x)-1]]),color=color))
            else:
                if i<n-1:
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[0],dots_up[i]],[dots_up[0],dots_low[1]],[dots_up[i+1],dots_low[0]]),color=color))
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[i],dots_up[0]],[dots_up[1],dots_low[0]],[dots_up[0],dots_low[i+1]]),color=color))
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[len(x)-1],dots_up[len(x)-1-i]],[dots_up[len(x)-1],dots_low[len(x)-2]],[dots_up[len(x)-2-i],dots_low[len(x)-1]]),color=color))
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_up[len(x)-1],dots_low[len(x)-1-i]],[dots_low[len(x)-1],dots_up[len(x)-2]],[dots_low[len(x)-2-i],dots_up[len(x)-1]]),color=color))
                else:
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[0],dots_up[i]],[dots_up[0],dots_low[1]],[dots_low[1],dots_low[0]]),color=color))
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[i],dots_up[0]],[dots_up[1],dots_low[0]],[dots_up[0],dots_up[1]]),color=color))
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[len(x)-1],dots_up[len(x)-1-i]],[dots_up[len(x)-1],dots_low[len(x)-2]],[dots_low[len(x)-2],dots_low[len(x)-1]]),color=color))
                    plt.gca().add_patch(plt.Polygon(intersect_3([dots_up[len(x)-1],dots_low[len(x)-1-i]],[dots_low[len(x)-1],dots_up[len(x)-2]],[dots_up[len(x)-2],dots_up[len(x)-1]]),color=color))

        for i in range(len(x)-3):
            m=Mach_from_omega(omega_n)
            color='#'+r+str(hex(int(k1*m+k0))).split('x')[-1]+b
            plt.gca().add_patch(plt.Polygon(intersect_3([dots_low[i+1],dots_low[i+2]],[dots_up[0],dots_low[i+2]],[dots_up[len(x)-1],dots_low[i+1]]),color=color))
            plt.gca().add_patch(plt.Polygon(intersect_3([dots_up[i+1],dots_up[i+2]],[dots_low[0],dots_up[i+2]],[dots_low[len(x)-1],dots_up[i+1]]),color=color))

        for i in range(n-1):
            for j in range(n-1-i):
                m=Mach_from_omega(omega(I(theta[i+1],omega_n),(j+2)*t_n/n+1000-I(theta[i],omega_n)))
                color='#'+r+str(hex(int(k1*m+k0))).split('x')[-1]+b
                plt.gca().add_patch(plt.Polygon(intersect_4([dots_up[0],dots_low[i+1]],[dots_low[0],dots_up[j+1+i]],[dots_up[0],dots_low[i+2]],[dots_low[0],dots_up[j+2+i]]),color=color))
                plt.gca().add_patch(plt.Polygon(intersect_4([dots_low[0],dots_up[i+1]],[dots_up[0],dots_low[j+1+i]],[dots_low[0],dots_up[i+2]],[dots_up[0],dots_low[j+2+i]]),color=color))
                plt.gca().add_patch(plt.Polygon(intersect_4([dots_up[len(x)-1],dots_low[len(x)-i-2]],[dots_low[len(x)-1],dots_up[len(x)-j-2-i]],[dots_up[len(x)-1],dots_low[len(x)-i-3]],[dots_low[len(x)-1],dots_up[len(x)-j-3-i]]),color=color))
                plt.gca().add_patch(plt.Polygon(intersect_4([dots_low[len(x)-1],dots_up[len(x)-i-2]],[dots_up[len(x)-1],dots_low[len(x)-j-2-i]],[dots_low[len(x)-1],dots_up[len(x)-i-3]],[dots_up[len(x)-1],dots_low[len(x)-j-3-i]]),color=color))
            m=Mach_from_omega(omega(I(theta[i+1],omega_n),t_n+1000-I(t_n,omega_n)))
            color='#'+r+str(hex(int(k1*m+k0))).split('x')[-1]+b
            plt.gca().add_patch(plt.Polygon(intersect_4([dots_up[0],dots_low[i+1]],[dots_low[0],dots_up[len(x)-2]],[dots_up[0],dots_low[i+2]],[dots_low[1],dots_up[len(x)-1]]),color=color))
            plt.gca().add_patch(plt.Polygon(intersect_4([dots_low[0],dots_up[i+1]],[dots_up[0],dots_low[len(x)-2]],[dots_low[0],dots_up[i+2]],[dots_up[1],dots_low[len(x)-1]]),color=color))
            plt.gca().add_patch(plt.Polygon(intersect_4([dots_up[len(x)-1],dots_low[len(x)-i-2]],[dots_low[len(x)-1],dots_up[1]],[dots_up[len(x)-1],dots_low[len(x)-i-3]],[dots_low[len(x)-2],dots_up[0]]),color=color))
            plt.gca().add_patch(plt.Polygon(intersect_4([dots_low[len(x)-1],dots_up[len(x)-i-2]],[dots_up[len(x)-1],dots_low[1]],[dots_low[len(x)-1],dots_up[len(x)-i-3]],[dots_up[len(x)-2],dots_low[0]]),color=color))
            for j in range(i):
                m=Mach_from_omega(omega(I(theta[i+1],omega_n),theta[1+j]+1000-I(theta[j+1],omega_n)))
                color='#'+r+str(hex(int(k1*m+k0))).split('x')[-1]+b
                plt.gca().add_patch(plt.Polygon(intersect_4([dots_up[0],dots_low[i+1]],[dots_low[1+j],dots_up[len(x)-1]],[dots_up[0],dots_low[i+2]],[dots_low[2+j],dots_up[len(x)-1]]),color=color))
                plt.gca().add_patch(plt.Polygon(intersect_4([dots_low[0],dots_up[i+1]],[dots_up[1+j],dots_low[len(x)-1]],[dots_low[0],dots_up[i+2]],[dots_up[2+j],dots_low[len(x)-1]]),color=color))
                

        m=Mach_from_omega(omega(I(-theta[0],omega_n),II(theta[0],omega_n))
        color='#'+r+str(hex(int(k1*m+k0))).split('x')[-1]+b
        plt.gca().add_patch(plt.Polygon(intersect_4([dots_low[0],dots_up[len(x)-2]],[dots_low[len(x)-1],dots_up[1]],[dots_up[len(x)-1],dots_low[1]],[dots_low[len(x)-2],dots_up[0]]),color=color))
        
    if check_mach_var.get()==1:
        distribution()
        
    #Shock waves---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def shock_lines():
        u=n
        if ratio > 1:
            for i in range(u):
                s_x=np.array([x[0],x[i+1]])
                s_y=np.array([y1[0],y2[i+1]])
                plt.plot(s_x,s_y,'r-',lw=1)
            for i in range(u):
                s_x=np.array([x[0],x[i+1]])
                s_y=np.array([y2[0],y1[i+1]])
                plt.plot(s_x,s_y,'r--',lw=1)
            for i in range(u):
                s_x=np.array([max(x),x[i+1]])
                s_y=np.array([max(y2),y1[i+1]])
                plt.plot(s_x,s_y,'b-',lw=1)
            for i in range(u):
                s_x=np.array([max(x),x[i+1]])
                s_y=np.array([min(y1),y2[i+1]])
                plt.plot(s_x,s_y,'b--',lw=1)
        else:
            for i in range(u):
                s_x=np.array([x[0],x[i+1]])
                s_y=np.array([y1[0],y2[i+1]])
                plt.plot(s_x,s_y,'b-',lw=1)
            for i in range(u):
                s_x=np.array([x[0],x[i+1]])
                s_y=np.array([y2[0],y1[i+1]])
                plt.plot(s_x,s_y,'b--',lw=1)
            for i in range(u):
                s_x=np.array([max(x),x[i+1]])
                s_y=np.array([min(y2),y1[i+1]])
                plt.plot(s_x,s_y,'r-',lw=1)
            for i in range(u):
                s_x=np.array([max(x),x[i+1]])
                s_y=np.array([max(y1),y2[i+1]])
                plt.plot(s_x,s_y,'r--',lw=1)
    if check_lines_var.get()==1:
        shock_lines()  

    #Wall-----------------------------------------------------------------------------------------------------------------------
    wall_x = np.array([-20,10,10])
    wall_y_under = np.array([0,0,-10])
    wall_y_above = np.array([80,80,90])
    plt.plot(wall_x,wall_y_under,'-k',lw=3)
    plt.plot(wall_x,wall_y_above,'-k',lw=3)

    plt.plot(0,0,'-k',lw=3,label='Solid boundary')
    plt.plot(0,0,'--k',lw=1,label='Slip boundary')
    plt.plot(0,0,'r-',lw=1,label='Left compression wave')
    plt.plot(0,0,'r--',lw=1,label='Right compression wave')
    plt.plot(0,0,'b-',lw=1,label='Left expansion wave')
    plt.plot(0,0,'b--',lw=1,label='Right expansion wave')

    if check_legend_var.get()==1:
        plt.legend(loc='lower left',labelspacing=0.1,fontsize='x-small',bbox_to_anchor=(0.05,0))

    title = 'Inlet M = '+str(int(M*100)/100)+', Outside/Inlet pressure ratio = '+str(int(ratio*100)/100)+', n = '+str(n)
    plt.title(title)

    plt.show()

window = tk.Tk()
window.title('Diamond')
window.geometry('220x260')
window.resizable(width=0, height=0)

entry_M = tk.Entry(width=20)
title_M = tk.Label(text="Inlet Mach number:")
title_M.pack()
entry_M.pack()

entry_ratio = tk.Entry(width=20)
title_ratio = tk.Label(text="Outside/inlet pressure ratio:")
title_ratio.pack()
entry_ratio.pack()

entry_n = tk.Spinbox(from_=2, to=10000,width=18)
title_n = tk.Label(text="Number of flow turns:")
title_n.pack()
entry_n.pack()

check_lines_var = tk.IntVar()
check_lines = tk.Checkbutton(text="Shock lines", variable=check_lines_var, onvalue=1, offvalue=0)
check_lines.pack()

check_mach_var = tk.IntVar()
check_mach = tk.Checkbutton(text="Mach numbers distribution", variable=check_mach_var, onvalue=1, offvalue=0)
check_mach.pack()

check_legend_var = tk.IntVar()
check_legend = tk.Checkbutton(text="Legend", variable=check_legend_var, onvalue=1, offvalue=0)
check_legend.pack()

start_button = tk.Button(
    text="Start visualisation",
    width=20,
    height=3, command=start
)
start_button.pack()

window.mainloop()
