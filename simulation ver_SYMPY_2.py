import math as ma
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import mpmath

#網格大小1*1*1km^3
#海、陸比熱
Co = 4200
Cl = 2500
#海、陸質量(暫時當做相等)(kg)
Mo=3e8
Ml=3e8

#模擬間隔
#t=int(input('interval(s) = '))
t=600
#模擬次數
times=round(86400*30/t)
clock=np.arange(1.,times,1.)
#短波輻射(W/m2)
S=1366
#輻射入射角度(度)
theta=0
#地表接收輻射
Sr = 0
#反射率
a = 0.3
#常數
e=5.67e-8
#畫圖用陣列
Toa=np.arange(1.,times,1.)
Tla=np.arange(1.,times,1.)
To=np.arange(1.,times,1.)
Tl=np.arange(1.,times,1.)
Top0=np.arange(1.,times,1.)
Top1=np.arange(1.,times,1.)
Ta2=np.arange(1.,times,1.)
dT=np.zeros(times-1)
Gh=np.zeros(times-1)
#初始溫度
T02=290#海洋第二層大氣
T12=290#陸地第二層大氣
T00=300.#海洋第一層大氣
T01=300.#陸地第一層大氣
T10=310.#海洋地表
T11=310.#陸地地表
#地表熱容量變化
dQo=0#海
dQl=0#陸
#對流公式:q=hA(Ta-Ts)
ho=0.6#水散熱係數
hl=2#岩石散熱係數
h=10#自然對流常數
#地表面積
Ao=1e6#海
Al=1e6#陸

#模擬程式
for i in range(times):
    #計算入射角度
    Sr=S*ma.sin(2*ma.pi*theta/86400)
    if Sr<0:
        Sr=0
        S0=0
    else:
        S0=1366
    #第二層大氣
    T02=round(((T00**4)/2)**0.25,4)#海
    T12=round(((T01**4)/2)**0.25,4)#陸

    dQo = round(Sr*(1-a)*Ao-Ao*e*T10**4+Ao*e*T00**4-ho*(T10-T00)*Ao,4)*t#計算水域熱量變化(以初始溫度時為零)
    T00=round(((T10**4)/2+ho*(T10-T00)/e/2-h*(T00-T01)/e/2+(T02**4)/2)**0.25,4)#以熱量平恆計算水域上方大氣溫度
    
    T10+= round(dQo/Mo/Co,4) #以熱量變化計算水域溫度變化
    
    dQl = round(Sr*(1-a)*Al-Al*e*T11**4+Al*e*T01**4-hl*(T11-T01)*Al,4)*t#計算陸地熱量變化(以初始溫度時為零)
    
    T01=round(((T11**4)/2+hl*(T11-T01)/e/2-h*(T01-T00)/e/2+(T12**4)/2)**0.25,4)#以熱量平恆計算陸地上方大氣溫度
    T11+= round(dQl/Ml/Cl,4)#以熱量變化計算陸地溫度變化

    theta+=1*t
    #將結果記錄在陣列中
    Toa[i-1]=T00
    Tla[i-1]=T01
    To[i-1]=T10
    Tl[i-1]=T11
    Ta2[i-1]=T12
    Top0[i-1]=T02
    dT[i-1]=T00-T01
    Gh[i-1]=-((8.314462618*dT[i-1])/(0.0289647*27750)) #計算氣壓梯度力(加速度) 假設兩地距離27.75km

#計算風速(m/s)
ws = np.zeros(times-1)
x = sp.symbols('x', commutative=True)
for i in range(times-1):
    x0, x1 = 0, i
    ws[i-1] = sp.integrate((Gh[i]), (x, x0, x1), manual=True).evalf()

plt.plot(clock,Toa,'g.')
plt.plot(clock,Tla,'b.')
plt.plot(clock,Gh,'y.')
plt.plot(clock,ws,'r.')
#plt.xlim([6000,times])
plt.show()