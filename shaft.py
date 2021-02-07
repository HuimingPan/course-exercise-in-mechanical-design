from math import ceil

class shaft():
    num=0
    def __init__(self,P,n):
        shaft.num+=1
        print("\n-----------------")
        print("Now, begin the design of shaft {}.\n".format(shaft.num))
        self.P=P
        self.n=n
        self.material='45'
        self.dia()
    def dia(self):
        if (shaft.num==1)&(self.material=='45'):
            self.A0=110
            self.tau_T_avaliable=28
            self.keyhole=1
        elif (shaft.num==2)&(self.material=='45'):
            self.A0=107
            self.tau_T_avaliable=34
            self.keyhole=1
        elif (shaft.num==3)&(self.material=='45'):
            self.A0=97
            self.tau_T_avaliable=40
            self.keyhole=1
        elif (shaft.num==4)&(self.material=='45'):
            self.A0=100
            self.tau_T_avaliable=40
            self.keyhole=1
        self.d0=self.A0*(self.P/self.n)**(1/3)
        self.d=ceil(self.d0*(1+0.1*self.keyhole))
        print('对于第 {} 根轴'.format(shaft.num))
        print('\t取A_0={},[τ_T]={}MPa'.format(self.A0,self.tau_T_avaliable))
        print('\t估算轴的直径 d={}mm'.format(self.d0))
        print('\t考虑键槽的影响，扩大轴的直径')
        print('\t    轴上有 {} 个键槽，直径增大{}%'.format(self.keyhole,self.keyhole*5))
        print('\t    则 d={}mm'.format(self.d))