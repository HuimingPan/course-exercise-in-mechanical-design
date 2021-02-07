from math import ceil
class rolling():
    num=0
    def __init__(self,P,n,d):
        rolling.num+=1
        print("\n-----------------")
        print("Now, begin the design of rolling bearing on shaft {}.\n".format(rolling.num))
        self.P=P
        self.n=n
        self.d_min=d
        self.T=9550*P/n
        print('{}.计算第{}根轴上的轴承'.format(rolling.num,rolling.num))
        print('\t1.滚动轴承的基本要求参数')
        print('\t  滚动轴承的输入功率 P={} kW'.format(self.P))
        print('\t  滚动轴承的转速 n={} r/min'.format(self.n))
        print('\t  滚动轴承公称转矩 T={} N·m'.format(self.T))
        print('\t  初算的轴的直径为 {} mm'.format(self.d_min))
        self.type_select()
        self.size_confirm()
    def type_select(self):
        self.type='圆锥滚子轴承'
        self.type_code='3'
        self.advantage='装拆调整方便、价格较低'
        print('\t2.滚动轴承的类型选择')
        print('\t  由于轴承受到轴向和径向力作用')
        print('\t  且{}具有{}的优点'.format(self.type,self.advantage))
        print('\t  故选用{}'.format(self.type))
    def size_confirm(self):
        self.inner_dia=ceil(self.d_min/5)*5
        self.inner_dia_code=int(self.inner_dia/5)
        self.external_dia_code=3
        self.width_code=0
        print('\t3.选择滚动轴承的尺寸')
        print('\t  滚动轴承的内径为{} mm,内径代号为{}'.format(self.inner_dia,self.inner_dia_code))
        print('\t  滚动轴承的直径代号为{}'.format(self.external_dia_code))
        print('\t  滚动轴承的宽度代号为{}'.format(self.width_code))
        print('\t4.最终确定轴承的代号为')
        print('\t  {}{}{}{:0>2}'.format(self.type_code,self.width_code,self.external_dia_code,self.inner_dia_code))