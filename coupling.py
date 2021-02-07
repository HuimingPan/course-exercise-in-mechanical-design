
class coupling():
    def __init__(self,P,n):
        print("\n-----------------")
        print("Now, begin the design of coupling.\n")
        self.P=P
        self.n=n
        self.T=9550*P/n
        print('联轴器的基本要求参数')
        print('\t联轴器的输入功率 P={} kW'.format(self.P))
        print('\t联轴器的转速 n={} r/min'.format(self.n))
        print('\t联轴器公称转矩 T={} N·m'.format(self.T))
        self.type_select()
        self.size_confirm()
    def type_select(self):
        self.type='弹性套柱销联轴器'
        print('选择联轴器类型')
        print('\t由于联轴器在减速器的输出端，转速低，传递转矩大，')
        print('\t故选用{}'.format(self.type))
    def size_confirm(self):
        self.K_A=1.5
        self.T_ca=self.K_A*self.T
        self.size='LT11'
        self.T_avalible=6300
        self.n_avaliable=1800
        self.d_min=80
        self.d_max=110
        print('确定联轴器型号')
        print('\t使用系数 K_A=',self.K_A)
        print('\t计算转矩 T_ca={} N·mm'.format(self.T_ca))
        print('\t选用{}型号弹性套柱销联轴器'.format(self.size))
        print('\t   许用转矩 [T]=',self.T_avalible,'N·m')
        print('\t   许用转速 [n]=',self.n_avaliable,'r/min')
        print('\t   轴颈在{}mm~{}mm之间'.format(self.d_min,self.d_max))
        print('\t故{}型号弹性套柱销联轴器适用'.format(self.size))