from math import pi,sin,ceil
class belt():
    '''
    '''
    def __init__(self,P,n,i):
        print("\n-----------------")
        print("Now,begin the design of belt\n")
        self.P=P
        self.n1=n
        self.i=i
        self.belt_design()
    
    def belt_design(self):
        print('1.确定计算功率 P_ca')
        self.K_A=1.3
        self.P_ca=self.K_A*self.P
        print('\t工作系数 K_A=',self.K_A)
        print('\t计算功率 P_ca={}kW'.format(self.P_ca))
        print('2.选择V带带型')
        self.type='A'
        print('\t根据P_ca、n1，选择{}型V带'.format(self.type))
        print('3.确定带轮的基准直径d_d,并验算带速 v')
        self.dd1=125#查表8-7、8-9、图8-11得
        self.v=pi*self.dd1*self.n1/60/1000
        self.dd20=self.dd1*self.i
        self.dd2=560#查表8-9 根据计算值取标准值
        print('\t初取小带轮基准直径 d_d1={}mm'.format(self.dd1))
        print('\t带速 v={}m/s'.format(self.v))
        if (self.v<30)&(self.v>5):
            print('\t带速在5m/s~30m/s之间，合适')
        else :
            print('\t带速不合适！！！')
        print('\t则大带轮基准直径 d_d2={}mm,取为 {}mm'.format(self.dd20,self.dd2))
        print('4.确定V带的中心距a和基准直径L_d')
        self.a0=0.7*(self.dd1+self.dd2)
        self.Ld0=2*self.a0+pi/2*(self.dd1+self.dd2)+(self.dd2-self.dd1)**2/4/self.a0
        self.Ld=2200
        self.a=self.a0+(self.Ld-self.Ld0)/2
        self.a_min=self.a-0.015*self.Ld
        self.a_max=self.a+0.030*self.Ld
        print('\t(1)初算中心距 a0={}mm'.format(self.a0))
        print('\t(2)初算带所需基准长度 L_d0={}mm,选为 {}mm'.format(self.Ld0,self.Ld))
        print('\t(3)实际中心距 a={:.6}mm，变化范围为{:.6}~{:.6}mm'.format(self.a,self.a_min,self.a_max))
        print('5.验算小带轮上的包角 α1')
        self.alpha_1=180-(self.dd2-self.dd1)*57.3/self.a
        print('\t小带轮上的包角 α1={:.6}°'.format(self.alpha_1))
        if self.alpha_1>120:
            print("\t小带轮上的包角 α1 合适！")
        print('6.计算带的根数 z')
        print('\t(1)计算单根V带的额定功率 Pr')
        self.P0=1.3816#插值法得到 1.3816kW
        self.Delta_P0=0.1516#插值法得到 1.1516kW
        self.K_alpha=0.86
        self.K_L=1.06#表8-2
        self.Pr=(self.P0+self.Delta_P0)*self.K_alpha*self.K_L
        print('\t   单根V带的基本额定功率 P0={}kW'.format(self.P0))  
        print('\t   单根V带的额定功率增量 ΔP0={}kW'.format(self.Delta_P0)) 
        print('\t   修正系数 K_α={},K_L={}'.format(self.K_alpha,self.K_L))
        print("\t   单根V带的额定功率为 {} kW".format(self.Pr))
        print('\t(2)计算V带的根数率 z')
        self.z0=self.P_ca/self.Pr
        self.z=ceil(self.z0)
        print("\t   带的根数率 z={},取{}根".format(self.z0,self.z))
        print("7.确定带的初拉力 F0")
        if self.type=='B':
            self.q=0.170
        elif self.type=='A':
            self.q=0.105
        self.F0=500*(2.5-self.K_alpha)*self.P_ca/(self.K_alpha*self.z*self.v)+self.q*self.v**2
        print('\t{}型带的单位长度质量q={}kg/m'.format(self.type,self.q))
        print("\t单根V带的初拉力 F0={} N".format(self.F0))
        print("8.确定带的压轴力 F_p")
        self.F_p=2*self.z*self.F0*sin(self.alpha_1/2*pi/180)
        print("\t带的压轴力 F_p={} N".format(self.F_p))

    def belt_wheel(self,shaft_d):
        小带轮=0
        大带轮=1
        self.wheel_dd=[self.dd1,self.dd2]
        self.wheel_type=[]
        for d_d in [self.dd1,self.dd2]:
            if self.type=='B':
                b_d=14.0
                h_amin=3.5
                h_fmin=10.8
                e=19
                f_min=11.5
                delta_min=7.5
                if d_d<=190:
                    phi=34
                else:
                    phi=38
                d_a=d_d+h_amin+0.5
                
            if d_d<=2.5*shaft_d:
                self.wheel_type.append('实心式')
            elif d_d<=300:
                self.wheel_type.append('腹板式')
            elif d_d>300:
                self.wheel_type.append('轮辐式')
        print('小带轮基准直径为 {} mm,采用 {} 结构'.format(self.wheel_dd[小带轮],self.wheel_type[小带轮]))
        print('小带轮基准直径为 {} mm,采用 {} 结构'.format(self.wheel_dd[大带轮],self.wheel_type[大带轮]))
