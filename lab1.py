#Лаба 1 по ВычМату. Вариант 2.
#f(x)=exp(x)*sin(x) [1;2]
import numpy as np
import math as m
import matplotlib.pyplot as mp

#--------------------------------------------------------------
#---functions---
def f(x0):
    return m.exp(x0) * m.sin(x0)
def pkra_f(x_p, h_p):
    return (f(x_p + h_p) - f(x_p)) / h_p
def lkra_f(x_l, h_l):
    return (f(x_l) - f(x_l - h_l)) / h_l
def ckra1_f(x_c1, h_c1):
    return (f(x_c1 + h_c1) - f(x_c1 - h_c1)) / (2*h_c1)
def ckra2_f(x_c2, h_c2):
    return (f(x_c2 - h_c2) - 2*f(x_c2) + f(x_c2 + h_c2)) / (h_c2**2)
#--------------------------------------------------------------
h = float(input("Введите шаг, h = "))
eps = float(input("Введите значение эпсилон, eps = "))
#-----------------------------------
x = [ii for ii in np.arange(1.0, 2.0, h)]
y = [f(x[iii]) for iii in range(len(x))]
#-----------------------------------------------------
def Pkra_Lkra(x_, h_):
    while True: #для пкра и лкра
        r_pkra = []
        r_lkra = []
        pkra = []
        lkra = []
   
        #пкра
        for i1 in range(0, len(x_)-1):
            pkra.append(pkra_f(x_[i1], h_))
        pkra.append(0)
    
        for i1 in range(0, len(x_)-1):
            r_pkra.append(abs((pkra_f(x_[i1], h_) - pkra_f(x_[i1], h_/2))/(2**1 - 1)))
        r_pkra.append(0)
        Rmax_pkra = max(r_pkra)
    
        #лкра
        lkra.append(0)
        for i2 in range(1, len(x_)):
            lkra.append(lkra_f(x_[i2], h_))
    
        r_lkra.append(0)
        for i2 in range(1, len(x_)):
            r_lkra.append(abs((lkra_f(x_[i2], h_) - lkra_f(x_[i2], h_/2))/(2**1 - 1)))
        Rmax_lkra = max(r_lkra)
    
        if eps >= Rmax_pkra and eps >= Rmax_lkra:
            return x_, h_, pkra, lkra, r_pkra, r_lkra, Rmax_pkra, Rmax_lkra
    
        x_.clear()
        pkra.clear()
        lkra.clear()
        r_pkra.clear()
        r_lkra.clear()
    
        h_ /= 2
        x_ = [i for i in np.arange(1.0, 2.0, h_)]
    
x_pkra_lkra, h1, pkra, lkra, r_pkra, r_lkra, Rmax_pkra, Rmax_lkra = Pkra_Lkra(x, h)
y_pkra_lkra = [f(x_pkra_lkra[i3]) for i3 in range(len(x_pkra_lkra))]
#-----------------------------------------------------
def Ckra_1(x_, h_):
    while True: #для цкра1
        r_ckra1 = []
        ckra1 = []
        leng = len(x_) - 1
        
        ckra1.append(0)
        for j in range(1, leng):
            ckra1.append(ckra1_f(x_[j], h_))
        ckra1.append(0)
            
        r_ckra1.append(0)
        for j in range(1, leng):
            r_ckra1.append(abs((ckra1_f(x_[j], h_) - ckra1_f(x_[j], h_/2))/(2**2 - 1)))
        r_ckra1.append(0)
        Rmax_ckra1 = max(r_ckra1)
                
        if (eps >= Rmax_ckra1):
            return x_, h_, ckra1, r_ckra1, Rmax_ckra1
                
        x_.clear()
        ckra1.clear()
        r_ckra1.clear()
        
        h_ /= 2
        x_ = [i4 for i4 in np.arange(1.0, 2.0, h_)]
        
x = [ii for ii in np.arange(1.0, 2.0, h)]
x_ckra1, h2, ckra1, r_ckra1, Rmax_ckra1 = Ckra_1(x, h)
#-----------------------------------------------------
def Ckra_2(x_, h_):
    while True:
        r_ckra2 = []
        ckra2 = []
        leng = len(x_) - 1
        
        ckra2.append(0)
        for i6 in range(1, leng):
            ckra2.append(ckra2_f(x_[i6], h_))
        ckra2.append(0)
            
        r_ckra2.append(0)
        for i6 in range(1, leng):
            r_ckra2.append(abs((ckra2_f(x_[i6], h_) - ckra2_f(x_[i6], h_/2))/(2**2 - 1)))
        r_ckra2.append(0)
        Rmax_ckra2 = max(r_ckra2)
                
        if (eps >= Rmax_ckra2):
            return x_, h_, ckra2, r_ckra2, Rmax_ckra2
                
        x_.clear()
        ckra2.clear()
        r_ckra2.clear()
                
        h_ /= 2
        x_ = [i7 for i7 in np.arange(1.0, 2.0, h_)]
        
x = [ii for ii in np.arange(1.0, 2.0, h)]
x_ckra2, h3, ckra2, r_ckra2, Rmax_ckra2 = Ckra_2(x, h)
#-----------------------------------------------------
M2 = 2.938
M3 = 19.589
M4 = 26.876
E = 2**(-26)
r1 = (M2 * h1) / 2.0 #для пкра
r2 = (2 * E) / h1
r = r1 + r2

h0_pkra_lkra = 2 * m.sqrt(E / M2)
h0_ckra1 = ((3 * E) / M3)**(1/3)
h0_ckra2 = 2 * ((3 * E) / M4)**(1/4)

#-----------работа с файлом------------------------------
file = open('output_table.txt','w')
try:
    #вывод пкра и лкравместе с х и у
    st = "{:^10.6}{:^10.6}{:^10.6}{:^10.6}\n".format("x", "y", "ПКРА", "ЛКРА")
    file.write(st)
    for i in range(len(x_pkra_lkra)):
        st = "{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n".format(x_pkra_lkra[i], y_pkra_lkra[i], pkra[i], lkra[i])
        file.write(st)
    file.write("\n\n")
   
    #вывод цкра1
    st = "{:^10.6}{:^10.6}\n".format("x", "ЦКРА1")
    file.write(st)
    for i in range(len(x_ckra1)):
        st = "{:10.6f}{:10.6f}\n".format(x_ckra1[i], ckra1[i])
        file.write(st)
    file.write("\n\n")
    
    #вывод цкра2
    st = "{:^10.6}{:^10.6}\n".format("x", "ЦКРА2")
    file.write(st)
    for i in range(len(x_ckra2)):
        st = "{:10.6f}{:10.6f}\n".format(x_ckra2[i], ckra2[i])
        file.write(st)
    
    st = "\n\nОптимальный шаг для ПКРА и ЛКРА h0 = " + str(h0_pkra_lkra) + "\n\n"
    file.write(st)
    st = "Оптимальный шаг для 1й ЦКРА h0 = " + str(h0_ckra1) + "\n\n"
    file.write(st)
    st = "Оптимальный шаг 2й ЦКРА h0 = " + str(h0_ckra2) + "\n\n"
    file.write(st)
    
    st = "Суммарная погрешность для ПКРА и ЛКРА r = " + str(r) + "\n\n"
    file.write(st)
    st = "Суммарная погрешность для ЦКРА1 r = " + str((M3*(h2**2))/6 + E/h2) + "\n\n"
    file.write(st)
    st = "Суммарная погрешность для ЦКРА2 r = " + str((M4*(h3**2))/12 + (4*E)/(h3**2)) + "\n\n"
    
    file.write(st)
    st = "Rmax для ПКРА R = " + str(Rmax_pkra) + "\n\n"
    file.write(st)
    st = "Rmax для ЛКРА R = " + str(Rmax_lkra) + "\n\n"
    file.write(st)
    st = "Rmax для 1й ЦКРА R = " + str(Rmax_ckra1) + "\n\n"
    file.write(st)
    st = "Rmax для 2й ЦКРА R = " + str(Rmax_ckra2) + "\n\n"
    file.write(st)
finally:
    file.close()

#---------------plot-------------------------------------
x = [ii for ii in np.arange(1.0, 2.0, h)]
y = [f(x[iii]) for iii in range(1, len(x)-1)]
y1 = [m.exp(x[i])*m.cos(x[i]) + m.exp(x[i])*m.sin(x[i]) for i in range(1, len(x)-1)]
y2 = [2*m.exp(x[i])*m.cos(x[i]) for i in range(1, len(x)-1)]
pkra = [pkra_f(x[i],h) for i in range(1, len(x)-1)]
lkra = [lkra_f(x[i],h) for i in range(1, len(x)-1)]
ckra1 = [ckra1_f(x[i],h) for i in range(1, len(x)-1)]
ckra2 = [ckra2_f(x[i],h) for i in range(1, len(x)-1)]

#обрезание концов у кра-шек
del x[0]
del x[len(x)-1]

mp.title("Построение графиков функций")
mp.xlabel("x")
mp.ylabel("y=f(x)")

mp.plot(x, y, color="red", label="y(x)", linestyle = ":", linewidth = 3)
mp.plot(x, y1, color="black", label="y\'(x)", linestyle = ":", linewidth = 6)
mp.plot(x, y2, color="blue", label="y\"(x)", linestyle = ":", linewidth = 6)
mp.plot(x, pkra, color="yellow", label="ПКРА", linestyle = ":", linewidth = 3)
mp.plot(x, lkra, color="orange", label="ЛКРА", linestyle = ":", linewidth = 3)
mp.plot(x, ckra1, color="green", label="ЦКРА1", linestyle = ":", linewidth = 3)
mp.plot(x, ckra2, color="pink", label="ЦКРА2", linestyle = ":", linewidth = 3)

mp.legend()
mp.show()
mp.savefig('output_pic.png')