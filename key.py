class key():
    num=0
    def __init__(self,T,d):
        self.T=T
        self.d=d
        key.num+=1
        self.key_type()
        self.key_size()
        self.load()
        print(self.T)
        print('σ_bs=',self.sigma_bs)
    def key_type(self):
        if (key.num==1)|(key.num==4):
            self.type='半圆普通平键'
        elif (key.num==2)|(key.num==3):
            self.type='矩形花键'

    def key_size(self):
        if (key.num==1)&(self.type=='半圆普通平键'):
            self.b=10
            self.h=8
            self.L=56
            self.l=self.L-0.5*self.b
        elif (key.num==2)&(self.type=='矩形花键'):
            self.N=8
            self.d=62
            self.D=72
            self.B=12
            
            self.C=0.6

            self.Psi=0.7
            self.l=self.B
            self.h=(self.D-self.d)/2-2*self.C
            self.d_m=(self.D+self.d)/2
        elif (key.num==3)&(self.type=='矩形花键'):
            self.N=10
            self.d=92
            self.D=102
            self.B=14
            
            self.C=0.4

            self.Psi=0.7
            self.l=self.B
            self.h=(self.D-self.d)/2-2*self.C
            self.d_m=(self.D+self.d)/2
        elif(key.num==4)&(self.type=='半圆普通平键'):
            self.b=25
            self.h=14
            self.L=280
            self.l=self.L-0.5*self.b

    def load(self):
        if self.type=='半圆普通平键':
            self.sigma_bs=4000*self.T/(self.h*self.l*self.d)
        elif self.type=='矩形花键':
            self.sigma_bs=2000*self.T/(self.Psi*self.N*self.h*self.l*self.d_m)
        