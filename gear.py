from math import pi,sin,cos,tan,atan,acos,sqrt,ceil,floor,radians,degrees

class gear():
    '''
    '''
    number=0
    def __init__(self,P,n,i):
        gear.number+=1
        print("\n-----------------")
        print("Now, begin the design of gear {}.\n".format(gear.number))
        self.P=P
        self.n1=n
        self.i=i
        self.T=9550*P/self.n1
        self.design(self.P,self.n1,self.i)
    def design(self,P,n,i):
        '''
        This function is programmed to design a gear according to “《机械设计》(第十版) 濮良贵等主编”
        We should first confirm these data：
            --input power P(kW)
            --rotate speed n(r/min)
            --input torque T
            --transmission rate i
            --the operating condition
        '''
        n1=n
        precision=7#variable precision represents the precision level.
        if gear.number==1:
            z1=23
            self.material={'小齿轮':{'牌号':'20Cr','热处理方法':'渗碳后淬火','齿面硬度':'58~62HRC','接触疲劳极限':1400,'弯曲疲劳极限':800},
                     '大齿轮':{'牌号':'20Cr','热处理方法':'渗碳后淬火','齿面硬度':'58~62HRC','接触疲劳极限':1400,'弯曲疲劳极限':800}
                   }
            #self.material={'小齿轮':{'牌号':'35SiMn','热处理方法':'调质','齿面硬度':'310~360HBW','接触疲劳极限':840,'弯曲疲劳极限':680},
             #         '大齿轮':{'牌号':'35SiMn','热处理方法':'调质','齿面硬度':'310~360HBW','接触疲劳极限':840,'弯曲疲劳极限':680}
            #        }
        else :
            z1=20
            self.material={'小齿轮':{'牌号':'20Cr2Ni4','热处理方法':'渗碳后淬火','齿面硬度':'58~62HRC','接触疲劳极限':1650,'弯曲疲劳极限':1050},
                      '大齿轮':{'牌号':'20Cr2Ni4','热处理方法':'渗碳后淬火','齿面硬度':'58~62HRC','接触疲劳极限':1650,'弯曲疲劳极限':1050}
                    }
        z2=ceil(z1*i)
        u=z2/z1
        self.alpha=20
        beta=14#variable beta represents the helix angel measured in degree.
        
        print("\n1.选择齿轮类型、精度等级、材料及齿数")
        print('\t(1)选用标准斜齿圆柱齿轮传动')
        print('\t(2)选用精度等级',precision)
        print('\t(3)材料选择:')
        print('\t   选择小齿轮材料为{}({})，齿面硬度{}'.format(self.material['小齿轮']['牌号'],self.material['小齿轮']['热处理方法'],self.material['小齿轮']['齿面硬度']))
        print('\t   选择大齿轮材料为{}({})，齿面硬度{}'.format(self.material['小齿轮']['牌号'],self.material['小齿轮']['热处理方法'],self.material['小齿轮']['齿面硬度']))
        print('\t(4)初选小齿轮齿数 ')
        print('\t   小齿轮、大齿轮齿数分别为 z1={},z2={}'.format(z1,z2))
        print('\t(5)初选螺旋角 β={}°，压力角 α={}°'.format(beta,self.alpha))

        print('\n2.按齿面接触疲劳强度设计')
        print('\t(1)由式(10-24)确定小齿轮分度圆直径')
        print('\t   1)确定公式中各参数值')
        K_Ht=1.3#试选载荷系数
        T1=9.55*(10**6)*P/n1
        self.alpha_t=atan(tan(radians(self.alpha))/cos(radians(beta)))
        beta_b=atan(tan(radians(beta))*cos(self.alpha_t))
        Z_H=sqrt(2*cos(beta_b)/cos(self.alpha_t)/sin(self.alpha_t))
        Z_E=189.8#查表得到材料的弹性影响系数
        h_an_star=1#齿高系数选1？？
        c_n_star=0.25
        if gear.number==1:
            phi_d=0.8#查表取齿宽系数
        else :
            phi_d=0.7#查表取齿宽系数
        alpha_at1=acos(z1*cos(self.alpha_t)/(z1+2*h_an_star*cos(radians(beta))))
        alpha_at2=acos(z2*cos(self.alpha_t)/(z2+2*h_an_star*cos(radians(beta))))
        epsilon_alpha=(z1*(tan(alpha_at1)-tan(self.alpha_t))+z2*(tan(alpha_at2)-tan(self.alpha_t)))/(2*pi)
        epsilon_beta=(phi_d*z1*tan(radians(beta)))/(pi)
        Z_epsilon=sqrt((4-epsilon_alpha)*(1-epsilon_beta)/3+epsilon_beta/epsilon_alpha)
        Z_beta=sqrt(cos(radians(beta)))
        sigma_Hlim1=self.material['小齿轮']['接触疲劳极限']
        sigma_Hlim2=self.material['大齿轮']['接触疲劳极限']#查表得到高速级大、小齿轮的接触疲劳极限
        j=1
        Lh=2*12*360*10
        N1=60*n1*j*Lh
        N2=N1/u
        if gear.number==1:
            K_HN1=0.92#查图10-19可得
            K_HN2=0.97#查图10-19可得
        else :
            K_HN1=0.97#查图10-19可得
            K_HN2=1.0#查图10-19可得
        failure_probability=0.01
        S=1
        simga_H_avilable_1=K_HN1*sigma_Hlim1/S
        simga_H_avilable_2=K_HN2*sigma_Hlim2/S
        sigma_H_avilable=min((simga_H_avilable_1,simga_H_avilable_2))
        print('\t     试选载荷系数 K_Ht=',K_Ht)
        print('\t     小齿轮传递的转矩 T1={} N·mm'.format(T1))
        print('\t     选择齿宽系数Φ_d=',phi_d)
        print('\t     计算区域系数 Z_H') 
        print('\t         α_t=',degrees(self.alpha_t))
        print('\t         β_b=',degrees(beta_b)) 
        print('\t         Z_H=',Z_H)
        print('\t     弹性影响系数 Z_E={} MPa^0.5'.format(Z_E))
        print('\t     计算重合度系数Z_ε')
        print('\t         α_at1=',degrees(alpha_at1))
        print('\t         α_at2=',degrees(alpha_at2))
        print('\t         ε_α=',epsilon_alpha)
        print('\t         ε_β=',epsilon_beta)
        print('\t         Z_ε=',Z_epsilon)
        print('\t     螺旋角系数Z_β=',Z_beta)
        print('\t     计算接触疲劳许用应力[σ_H]')
        print('\t         查图得接触疲劳极限为 σ_Hlim1={} MPa,σ_Hlim2={} MPa'.format(sigma_Hlim1,sigma_Hlim2))
        print('\t         应力循环次数 N1={},N2={}'.format(N1,N2))
        print('\t         接触疲劳寿命系数 K_HN1={},K_HN2={}'.format(K_HN1,K_HN2))
        print('\t         取失效概率 {}%,安全系数 S={}'.format(failure_probability*100,S))
        print('\t         [σ_H1]={:.8} MPa,[σ_H2]={:.8} MPa'.format(simga_H_avilable_1,simga_H_avilable_2))
        print('\t         [σ_H]=min([σ_H1],[σ_H2])={:.8} MPa'.format(sigma_H_avilable))

        print('\t   2)试算小齿轮分度圆直径')
        d_1t_min=((2*K_Ht*T1/phi_d)*(1+1/u)*(Z_H*Z_E*Z_epsilon*Z_beta/sigma_H_avilable)**2)**(1/3)
        print('\t     d_1t>={:.8}mm'.format(d_1t_min))

        print('\t(2)调整小齿轮分度圆直径')
        print('\t   1)计算实际载荷系数前的数据准备')
        v=pi*d_1t_min*n1/60/1000
        b=phi_d*d_1t_min
        print('\t     圆周速度v={:.8}m/s'.format(v))
        print('\t     齿宽b={:.8}mm'.format(b))
        print('\t   2)计算实际载荷系数K_H')
        K_A=1.25#使用系数，查表10-2取
        K_V=1#查图10-8取
        F_t1=2*T1/d_1t_min
        if K_A*F_t1/b>=100:
            K_H_alpha=1.2#对于7级精度的选择
        else:
            K_H_alpha=1.4#
        K_H_beta=1.289#查表10-4，插值法得到0.003+1.287+(30.267-40)*(1.287-1.293)/(40-80)
        K_H=K_A*K_V*K_H_alpha*K_H_beta
        d_1H=d_1t_min*(K_H/K_Ht)**(1/3)
        m_nH=d_1H*cos(radians(beta))/z1
        print('\t     使用系数 KA=',K_A)
        print('\t     动载系数 KV=',K_V)
        print('\t     计算齿间载荷分布系数 K_Hα')
        print('\t         小齿轮的径向力 F_t1={:.8} N'.format(F_t1))
        print('\t         K_A*F_t1/b={:.8} N'.format(K_A*F_t1/b))
        print('\t         齿间载荷分布系数 K_Hα=',K_H_alpha)
        print('\t         齿向载荷分布系数 K_Hβ=',K_H_beta)
        print('\t     实际载荷系数 K_H=',K_H)
        print('\t   3)使用实际载荷系算得的分度圆直径 d_1H={:.5}mm'.format(d_1H))
        print('\t     相应的齿轮模数 m_nH={:.5}mm'.format(m_nH))

        print('\n3.按齿面弯曲疲劳强度设计')
        print('\t(1)按式(10-17)试算模数')
        print('\t   1)确定公式中的各种参数')
        K_Ft=1.3
        epsilon_alpha_v=epsilon_alpha/(cos(beta_b))**2
        Y_epsilon=0.25+0.75/epsilon_alpha_v
        Y_beta=1-epsilon_beta*beta/120
        z_v1=z1/(cos(radians(beta)))**3
        z_v2=z2/(cos(radians(beta)))**3
        if gear.number==1:
            Y_Fa1=2.62#查表10-5得到
            Y_Fa2=2.13#查表10-5得到
            Y_Sa1=1.59#查表10-5得到
            Y_Sa2=1.82#查表10-5得到
            K_FN1=0.90#查图10-18得
            K_FN2=0.95#查图10-18得
        else :
            Y_Fa1=2.72#查表10-5得到
            Y_Fa2=2.21#查表10-5得到
            Y_Sa1=1.57#查表10-5得到
            Y_Sa2=1.78#查表10-5得到
            K_FN1=0.98#查图10-18得
            K_FN2=0.99#查图10-18得
        sigma_Flim1=self.material['小齿轮']['弯曲疲劳极限']#查图10-20c得
        sigma_Flim2=self.material['大齿轮']['弯曲疲劳极限']#查图10-20c得
        S=1.4#弯曲疲劳安全系数
        sigma_F_avaliable_1=K_FN1*sigma_Flim1/S
        sigma_F_avaliable_2=K_FN2*sigma_Flim2/S
        num1=Y_Fa1*Y_Sa1/sigma_F_avaliable_1
        num2=Y_Fa2*Y_Sa2/sigma_F_avaliable_2
        num=max(num1,num2)
        print('\t     试选载荷系数 K_Ft=',K_Ft)
        print('\t     计算重合度系数 Y_ε') 
        print('\t         ε_αv=',epsilon_alpha_v)
        print('\t         重合度系数 Y_ε=',Y_epsilon) 
        print('\t     螺旋角系数 Y_β=',Y_beta)
        print('\t     计算(Y_Fa*Y_Sa/[σ_F])')
        print('\t         当量齿数 z_v1={:.8},z_v2={:.8}'.format(z_v1,z_v2))
        print('\t         齿形修正系数 Y_Fa1={:.8},Y_Fa2={:.8}'.format(Y_Fa1,Y_Fa2)) 
        print('\t         应力修正系数 Y_Sa1={:.8},Y_Sa2={:.8}'.format(Y_Sa1,Y_Sa2)) 
        print('\t         弯曲疲劳寿命系数 K_FN1={:.8},K_FN2={:.8}'.format(Y_Sa1,Y_Sa2))
        print('\t         弯曲疲劳安全系数 S=',S)
        print('\t         弯曲应力许用值 [σ_F]_1={:.8} MPa,[σ_F]_2={:.8} MPa'.format(sigma_F_avaliable_1,sigma_F_avaliable_2))
        print('\t         Y_Fa1*Y_Sa1/[σ_F1]={:.5},Y_Fa2*Y_Sa2/[σ_F2]={:.5}'.format(num1,num2))
        print('\t         Y_Fa*Y_Sa/[σ_F]=max(Y_Fa1*Y_Sa1/[σ_F1],Y_Fa2*Y_Sa2/[σ_F2])={:.5}'.format(num))
        print('\t   2)试算齿轮模数')
        m_nt_min=((2*K_Ft*T1*Y_epsilon*Y_beta*(cos(radians(beta)))**2)/(phi_d*z1**2)*num)**(1/3)
        print('\t     m_nt_min=',m_nt_min)
        print('\t(2)调整齿轮模数')
        print('\t   1)计算实际载荷系数前的数据准备')
        d1=m_nt_min*z1/cos(radians(beta))
        v=pi*d1*n1/60/1000
        b=phi_d*d1
        h=(2*h_an_star+c_n_star)*m_nt_min
        aspect_ratio=b/h
        print('\t     计算圆周速度 v')
        print('\t     圆周速度v={:.3}m/s'.format(v))
        print('\t     齿宽b={:.5}mm'.format(b))
        print('\t     宽高比 b/h')
        print('\t         高 h=',h)
        print('\t         宽高比 b/h=',aspect_ratio)
        print('\t   2)计算实际载荷系数')
        K_A=1
        K_V=1
        F_t1=2*T1/d1
        if K_A*F_t1/b>100:
            K_F_alpha=1.2
        else:
            K_F_alpha=1.4
        if gear.number==1:
            K_H_beta=1.289#查表10-4得
            K_F_beta=1.28#查图 10-13得
        else :
            K_H_beta=1.290##查表10-4，插值法得到0.003+1.287+(30.267-40)*(1.287-1.293)/(40-80)
            K_F_beta=1.22#查图 10-13得
        K_F=K_A*K_V*K_F_alpha*K_F_beta
        print('\t     动载系数 K_A=',K_A)
        print('\t     计算齿间载荷分布系数 K_Fα')
        print('\t         F_t1=',F_t1)
        print('\t         K_A*F_t1/b =',K_A*F_t1/b)
        print('\t         齿间载荷分布系数 K_Fα=',K_F_alpha)
        print('\t     齿间载荷分布系数 K_Hβ={},K_Fβ={}'.format(K_H_beta,K_F_beta))
        print('\t     实际载荷系数 K_F=',K_F)
        print('\t   3)按实际载荷系数算得齿轮模数和小齿轮分度圆直径')
        m_nF=m_nt_min*(K_F/K_Ft)**(1/3)
        d_1F=m_nF*z1/cos(radians(beta))
        if gear.number==1:
            self.m_n=2.0
        else:
            self.m_n=3.0
        d1=d_1H
        z10=d1*cos(radians(beta))/self.m_n
        z1=ceil(z10)
        z2=floor(z1*u)
        if gear.number==1:
            z2=floor(z1*u)
        else:
            z2=ceil(z1*u)
        print('\t     齿轮模数 m_nF=',m_nF)
        print('\t     小齿轮分度圆直径 d_1F=',d_1F)
        print('\t     对比按弯曲疲劳强度和接触疲劳强度的计算结果')
        print('\t         从满足齿根弯曲疲劳强度出发，就近圆整齿轮模数 m_n={}mm'.format(self.m_n))
        print('\t         按满足齿面解除疲劳强度出发，算得到分度圆直径 d1={:.5}mm'.format(d1))
        print('\t     则小齿轮齿数 z1={}，取为 z1={}'.format(z10,z1))
        print('\t     则大齿轮齿数 z2={}，取为 z2={}'.format(z1*u,z2))

        print('\n4.几何尺寸计算')
        print('\t(1)计算中心距')
        a0=(z1+z2)*self.m_n/2/cos(radians(beta))
        self.a=ceil(a0)
        if gear.number==2:
            self.a=self.a-1
        print('\t      计算得中心距 a=',a0)
        print('\t      圆整得中心距 a=',self.a)
        print('\t(2)按圆整后的中心距修正螺旋角')
        self.beta=degrees(acos((z1+z2)*self.m_n/2/self.a))
        print('\t      修正后的螺旋角 β={}°'.format(self.beta))
        print('\t(3)计算大小齿轮的分度圆直径')
        self.d1=z1*self.m_n/cos(radians(beta))
        self.d2=z2*self.m_n/cos(radians(beta))
        self.d_a1=self.d1+2*self.m_n
        self.d_a2=self.d2+2*self.m_n
        print('\t      小齿轮分度圆直径 d1={}'.format(self.d1))
        print('\t      大齿轮分度圆直径 d2={}'.format(self.d2))
        print('\t      小齿轮分度圆直径 d1={}'.format(self.d1))
        print('\t      大齿轮分度圆直径 d2={}'.format(self.d2))
        print('\t(4)计算齿轮宽度')
        b=phi_d*d1
        self.b1=ceil(b+5)
        self.b2=ceil(b)
        print('\t      设计齿轮宽度 b={}mm'.format(b))
        print('\t      取小齿轮宽度 b1={}mm'.format(self.b1))
        print('\t      取大齿轮宽度 b2={}mm'.format(self.b2))


# def gear_design(P,n,i):
#     '''
#     This function is programmed to design a gear according to “《机械设计》(第十版) 濮良贵等主编”
#     We should first confirm these data：
#         --input power P(kW)
#         --rotate speed n(r/min)
#         --input torque T
#         --transmission rate i
#         --the operating condition
#     '''
#     print("\n-----------------")
#     print("Now, begin the design of gear.\n")
#     n1=n
#     precision=7#variable precision represents the precision level.
#     z1=23
#     z2=ceil(z1*i)
#     u=z2/z1
#     alpha=20
#     beta=14#variable beta represents the helix angel measured in degree.
#     print("\n1.选择齿轮类型、精度等级、材料及齿数")
#     print('\t(1)选用标准斜齿圆柱齿轮传动')
#     print('\t(2)选用精度等级',precision)
#     print('\t(3)材料选择，选择小齿轮材料为40Cr(调质)，齿面硬度280HBW；\n\t   大齿轮材料选用45钢(调质)，齿面硬度240HBW')
#     print('\t(4)初选小齿轮齿数 ')
#     print('\t   小齿轮、大齿轮齿数分别为 z1={},z2={}'.format(z1,z2))
#     print('\t(5)初选螺旋角 β=14°，压力角 α=20°')

#     print('\n2.按齿面接触疲劳强度设计')
#     print('\t(1)由式(10-24)确定小齿轮分度圆直径')
#     print('\t   1)确定公式中各参数值')
#     K_Ht=1.3#试选载荷系数
#     T1=9.55*(10**6)*P/n1
#     alpha_t=atan(tan(radians(alpha))/cos(radians(beta)))
#     beta_b=atan(tan(radians(beta))*cos(alpha_t))
#     Z_H=sqrt(2*cos(beta_b)/cos(alpha_t)/sin(alpha_t))
#     Z_E=189.8#查表得到材料的弹性影响系数
#     h_an_star=1#齿高系数选1？？
#     c_n_star=0.25
#     phi_d=0.8#查表取齿宽系数
#     alpha_at1=acos(z1*cos(alpha_t)/(z1+2*h_an_star*cos(radians(beta))))
#     alpha_at2=acos(z2*cos(alpha_t)/(z2+2*h_an_star*cos(radians(beta))))
#     epsilon_alpha=(z1*(tan(alpha_at1)-tan(alpha_t))+z2*(tan(alpha_at2)-tan(alpha_t)))/(2*pi)
#     epsilon_beta=(phi_d*z1*tan(radians(beta)))/(pi)
#     Z_epsilon=sqrt((4-epsilon_alpha)*(1-epsilon_beta)/3+epsilon_beta/epsilon_alpha)
#     Z_beta=sqrt(cos(radians(beta)))
#     sigma_Hlim1=600
#     sigma_Hlim2=550#查表得到大小齿轮的接触疲劳极限 
#     j=1
#     Lh=2*12*360*10
#     N1=60*n1*j*Lh
#     N2=N1/u
#     K_HN1=0.92#查图10-19可得
#     K_HN2=0.97#查图10-19可得
#     failure_probability=0.01
#     S=1
#     simga_H_avilable_1=K_HN1*sigma_Hlim1/S
#     simga_H_avilable_2=K_HN2*sigma_Hlim2/S
#     sigma_H_avilable=min((simga_H_avilable_1,simga_H_avilable_2))
#     print('\t     试选载荷系数 K_Ht=',K_Ht)
#     print('\t     小齿轮传递的转矩T1=',T1)
#     print('\t     选择齿宽系数Φ_d=',phi_d)
#     print('\t     计算区域系数 Z_H') 
#     print('\t         α_t=',degrees(alpha_t))
#     print('\t         β_b=',degrees(beta_b)) 
#     print('\t         Z_H=',Z_H)
#     print('\t     弹性影响系数 Z_E=',Z_E)
#     print('\t     计算重合度系数Z_ε')
#     print('\t         α_at1=',degrees(alpha_at1))
#     print('\t         α_at2=',degrees(alpha_at2))
#     print('\t         ε_α=',epsilon_alpha)
#     print('\t         ε_β=',epsilon_beta)
#     print('\t         Z_ε=',Z_epsilon)
#     print('\t     螺旋角系数Z_β=',Z_beta)
#     print('\t     计算接触疲劳许用应力[σ_H]=',sigma_H_avilable,'MPa')
#     print('\t         查图得接触疲劳极限为 σ_Hlim1={} MPa,σ_Hlim2={} MPa'.format(sigma_Hlim1,sigma_Hlim2))
#     print('\t         应力循环次数 N1={},N2={}'.format(N1,N2))
#     print('\t         接触疲劳寿命系数 K_HN1={},K_HN2={}'.format(K_HN1,K_HN2))
#     print('\t         取失效概率 {}%,安全系数 S={}'.format(failure_probability*100,S))
#     print('\t         [σ_H1]={} MPa,[σ_H2]={} MPa'.format(simga_H_avilable_1,simga_H_avilable_2))
#     print('\t         [σ_H]=min([σ_H1],[σ_H2])={} MPa'.format(sigma_H_avilable))

#     print('\t   2)试算小齿轮分度圆直径')
#     d_1t_min=((2*K_Ht*T1/phi_d)*(1+1/u)*(Z_H*Z_E*Z_epsilon*Z_beta/sigma_H_avilable)**2)**(1/3)
#     print('\t     d_1t>={}mm'.format(d_1t_min))

#     print('\t(2)调整小齿轮分度圆直径')
#     print('\t   1)计算实际载荷系数前的数据准备')
#     v=pi*d_1t_min*n1/60/1000
#     b=phi_d*d_1t_min
#     print('\t     圆周速度v={}mm/s'.format(v))
#     print('\t     齿宽b={}mm'.format(v))
#     print('\t   2)计算实际载荷系数K_H')
#     K_A=1.25#使用系数，查表10-2取
#     K_V=1#查图10-8取
#     F_t1=2*T1/d_1t_min
#     if K_A*F_t1/b>=100:
#         K_H_alpha=1.2#对于7级精度的选择
#     else:
#         K_H_alpha=1.4#
#     K_H_beta=1.320#查表10-34得到
#     K_H=K_A*K_V*K_H_alpha*K_H_beta
#     d_1H=d_1t_min*(K_H/K_Ht)**(1/3)
#     m_nH=d_1H*cos(radians(beta))/z1
#     print('\t     使用系数 KA=',K_A)
#     print('\t     动载系数 KV=',K_V)
#     print('\t     计算齿间载荷分布系数 K_Hα=')
#     print('\t         小齿轮的径向力 F_t1=',F_t1)
#     print('\t         齿间载荷分布系数 K_Hα=',K_H_alpha)
#     print('\t     实际使用系数 K_H=',K_H)
#     print('\t     使用实际载荷系算得的分度圆直径 d_1H={}mm'.format(d_1H))
#     print('\t     相应的齿轮模数 m_nH=',m_nH)

#     print('\n3.按齿面弯曲疲劳强度设计')
#     print('\t(1)按式(10-17)试算模数')
#     print('\t   1)确定公式中的各种参数')
#     K_Ft=1.3
#     epsilon_alpha_v=epsilon_alpha/(cos(beta_b))**2
#     Y_epsilon=0.25+0.75/epsilon_alpha_v
#     Y_beta=1-epsilon_beta*beta/120
#     z_v1=z1/(cos(radians(beta)))**3
#     z_v2=z2/(cos(radians(beta)))**3
#     Y_Fa1=2.65#查表10-5得到
#     Y_Fa2=2.13#查表10-5得到
#     Y_Sa1=1.59#查表10-5得到
#     Y_Sa2=1.82#查表10-5得到
#     sigma_Flim1=500#查图10-20c得
#     sigma_Flim2=320#查图10-20c得
#     K_FN1=0.90#查图10-18得
#     K_FN2=0.95#查图10-18得
#     S=1.4#弯曲疲劳安全系数
#     sigma_F_avaliable_1=K_FN1*sigma_Flim1/S
#     sigma_F_avaliable_2=K_FN2*sigma_Flim2/S
#     num1=Y_Fa1*Y_Sa1/sigma_F_avaliable_1
#     num2=Y_Fa2*Y_Sa2/sigma_F_avaliable_2
#     num=max(num1,num2)
#     print('\t     试选载荷系数 K_Ft=',K_Ft)
#     print('\t     计算重合度系数 Y_ε') 
#     print('\t         ε_αv=',epsilon_alpha_v)
#     print('\t         重合度系数 Y_ε=',Y_epsilon) 
#     print('\t     螺旋角系数 Y_β=',Y_beta)
#     print('\t     计算(Y_Fa*Y_Sa/[σ_F])')
#     print('\t         当量齿数 Z_v1={},Z_v2={}'.format(z_v1,z_v2))
#     print('\t         齿形修正系数 Y_Fa1={},Y_Fa2={}'.format(Y_Fa1,Y_Fa2)) 
#     print('\t         应力修正系数 Y_Sa1={},Y_Sa2={}'.format(Y_Sa1,Y_Sa2)) 
#     print('\t         弯曲疲劳寿命系数 K_FN1={},K_FN2={}'.format(Y_Sa1,Y_Sa2))
#     print('\t         弯曲疲劳安全系数 S=',S)
#     print('\t         弯曲应力许用值 [σ_F]_1={},[σ_F]_2={}'.format(sigma_F_avaliable_1,sigma_F_avaliable_2))
#     print('\t         Y_Fa1*Y_Sa1/[σ_F1]={},Y_Fa2*Y_Sa2/[σ_F2]={}'.format(num1,num2))
#     print('\t         Y_Fa*Y_Sa/[σ_F]=max(Y_Fa1*Y_Sa1/[σ_F1],Y_Fa2*Y_Sa2/[σ_F2])={}'.format(num))
#     print('\t   2)试算齿轮模数')
#     m_nt_min=((2*K_Ft*T1*Y_epsilon*Y_beta*(cos(radians(beta)))**2)/(phi_d*z1**2)*num)**(1/3)
#     print('\t     m_nt_min=',m_nt_min)
#     print('\t(2)调整齿轮模数')
#     print('\t   1)计算实际载荷系数前的数据准备')
#     d1=m_nt_min*z1/cos(radians(beta))
#     v=pi*d1*n1/60/1000
#     b=phi_d*d1
#     h=(2*h_an_star+c_n_star)*m_nt_min
#     aspect_ratio=b/h
#     print('\t     计算圆周速度 v')
#     print('\t         小齿轮分度圆直径 d1=',d1)
#     print('\t         圆周速度 v=',v) 
#     print('\t     齿宽 b=',b)
#     print('\t     宽高比 b/h')
#     print('\t         高 h=',h)
#     print('\t         宽高比 b/h=',aspect_ratio)
#     print('\t   2)计算实际载荷系数')
#     K_A=1
#     K_V=1
#     F_t1=2*T1/d1
#     if K_A*F_t1/b>100:
#         K_F_alpha=1.2
#     else:
#         K_F_alpha=1.4
#     K_H_beta=1.200#查表10-4得
#     K_F_beta=1.19#查图 10-13得
#     K_F=K_A*K_V*K_F_alpha*K_F_beta
#     print('\t     动载系数 K_A',K_A)
#     print('\t     计算齿间载荷分布系数 K_Fα')
#     print('\t         F_t1=',F_t1)
#     print('\t         K_A*F_t1/b =',K_A*F_t1/b)
#     print('\t         齿间载荷分布系数 K_Fα=',K_F_alpha)
#     print('\t     齿间载荷分布系数 K_Hβ={},K_Fβ={}'.format(K_H_beta,K_F_beta))
#     print('\t     实际载荷系数 K_F=',K_F)
#     print('\t   3)按实际载荷系数算得齿轮模数和小齿轮分度圆直径')
#     m_nF=m_nt_min*(K_F/K_Ft)**(1/3)
#     d_1F=m_nF*z1/cos(radians(beta))
#     m_n=ceil(max(m_nF,m_nH))
#     d1=min(d_1F,d_1H)
#     z10=d_1H*cos(radians(beta))/m_n
#     z1=ceil(z1)
#     z2=floor(z1*u)
#     print('\t     齿轮模数 m_nF=',m_nF)
#     print('\t     小齿轮分度圆直径 d_1F=',d_1F)
#     print('\t     对比按弯曲疲劳强度和接触疲劳强度的计算结果')
#     print('\t         圆整齿轮模数 m_n=',m_n)
#     print('\t         分度圆直径 d1=',d1)
#     print('\t     则小齿轮齿数 z1={}，取为 z1={}'.format(z10,z1))
#     print('\t     则大齿轮齿数 z2={}，取为 z2={}'.format(z1*u,z2))

#     print('\n4.几何尺寸计算')
#     print('\t(1)计算中心距')
#     a0=(z1+z2)*m_n/2/cos(radians(beta))
#     a=ceil(a0)
#     print('\t      计算得中心距 a=',a0)
#     print('\t      圆整得中心距 a=',a)
#     print('\t(2)按圆整后的中心距修正螺旋角')
#     beta=degrees(acos((z1+z2)*m_n/2/a))
#     print('\t      修正后的螺旋角 β={}°'.format(beta))
#     print('\t(3)计算大小齿轮的分度圆直径')
#     d1=z1*m_n/cos(radians(beta))
#     d2=z2*m_n/cos(radians(beta))
#     print('\t      小齿轮分度圆直径 d1={}'.format(d1))
#     print('\t      大齿轮分度圆直径 d2={}'.format(d2))
#     print('\t(4)计算齿轮宽度')
#     b=phi_d*d1
#     b1=ceil(b+5)
#     b2=ceil(b)
#     print('\t      设计齿轮宽度 b={}mm'.format(b))
#     print('\t      取小齿轮宽度 b1={}mm'.format(b1))
#     print('\t      取大齿轮宽度 b2={}mm'.format(b2))

# def gear_design2(P,n,i):
#     '''
#     This function is programmed to design a gear according to “《机械设计》(第十版) 濮良贵等主编”
#     We should first confirm these data：
#         --input power P(kW)
#         --rotate speed n(r/min)
#         --input torque T
#         --transmission rate i
#         --the operating condition
#     '''
#     print("\n-----------------")
#     print("Now, begin the design of gear2.\n")
#     n1=n
#     precision=7#variable precision represents the precision level.
#     z1=20
#     z2=ceil(z1*i)
#     u=z2/z1
#     alpha=20#variable alpha represents the pressure angel measured in degree.
#     beta=14#variable beta represents the helix angel measured in degree.
#     print("\n1.选择齿轮类型、精度等级、材料及齿数")
#     print('\t(1)选用标准斜齿圆柱齿轮传动')
#     print('\t(2)选用精度等级',precision)
#     print('\t(3)材料选择，选择小齿轮材料为40Cr(调质)，齿面硬度280HBW；\n\t   大齿轮材料选用45钢(调质)，齿面硬度240HBW')
#     print('\t(4)初选小齿轮齿数 ')
#     print('\t   小齿轮、大齿轮齿数分别为 z1={},z2={}'.format(z1,z2))
#     print('\t(5)初选螺旋角 β=14°，压力角 α=20°')

#     print('\n2.按齿面接触疲劳强度设计')
#     print('\t(1)由式(10-24)确定小齿轮分度圆直径')
#     print('\t   1)确定公式中各参数值')
#     K_Ht=1.3#试选载荷系数
#     T1=9.55*(10**6)*P/n1
#     alpha_t=atan(tan(radians(alpha))/cos(radians(beta)))
#     beta_b=atan(tan(radians(beta))*cos(alpha_t))
#     Z_H=sqrt(2*cos(beta_b)/cos(alpha_t)/sin(alpha_t))
#     Z_E=189.8#查表得到材料的弹性影响系数
#     h_an_star=1#齿高系数选1？？
#     c_n_star=0.25
#     phi_d=0.7#查表取齿宽系数
#     alpha_at1=acos(z1*cos(alpha_t)/(z1+2*h_an_star*cos(radians(beta))))
#     alpha_at2=acos(z2*cos(alpha_t)/(z2+2*h_an_star*cos(radians(beta))))
#     epsilon_alpha=(z1*(tan(alpha_at1)-tan(alpha_t))+z2*(tan(alpha_at2)-tan(alpha_t)))/(2*pi)
#     epsilon_beta=(phi_d*z1*tan(radians(beta)))/(pi)
#     Z_epsilon=sqrt((4-epsilon_alpha)*(1-epsilon_beta)/3+epsilon_beta/epsilon_alpha)
#     Z_beta=sqrt(cos(radians(beta)))
#     sigma_Hlim1=1650
#     sigma_Hlim2=1650#查表得到大小齿轮的接触疲劳极限 
#     j=1
#     Lh=2*12*360*10
#     N1=60*n1*j*Lh
#     N2=N1/u
#     K_HN1=1.14#查图10-19可得
#     K_HN2=1.27#查图10-19可得
#     failure_probability=0.01
#     S=1
#     simga_H_avilable_1=K_HN1*sigma_Hlim1/S
#     simga_H_avilable_2=K_HN2*sigma_Hlim2/S
#     sigma_H_avilable=min((simga_H_avilable_1,simga_H_avilable_2))
#     print('\t    试选载荷系数 K_Ht=',K_Ht)
#     print('\t    小齿轮传递的转矩T1=',T1)
#     print('\t    选择齿宽系数Φ_d=',phi_d)
#     print('\t    计算区域系数 Z_H') 
#     print('\t        α_t=',degrees(alpha_t))
#     print('\t        β_b=',degrees(beta_b)) 
#     print('\t        Z_H=',Z_H)
#     print('\t    弹性影响系数 Z_E=',Z_E)
#     print('\t    计算重合度系数Z_ε')
#     print('\t        α_at1=',degrees(alpha_at1))
#     print('\t        α_at2=',degrees(alpha_at2))
#     print('\t        ε_α=',epsilon_alpha)
#     print('\t        ε_β=',epsilon_beta)
#     print('\t        Z_ε=',Z_epsilon)
#     print('\t    螺旋角系数Z_β=',Z_beta)
#     print('\t    计算接触疲劳许用应力[σ_H]=',sigma_H_avilable,'MPa')
#     print('\t        查图得接触疲劳极限为 σ_Hlim1={} MPa,σ_Hlim2={} MPa'.format(sigma_Hlim1,sigma_Hlim2))
#     print('\t        应力循环次数 N1={},N2={}'.format(N1,N2))
#     print('\t        接触疲劳寿命系数 K_HN1={},K_HN2={}'.format(K_HN1,K_HN2))
#     print('\t        取失效概率 {}%,安全系数 S={}'.format(failure_probability*100,S))
#     print('\t        [σ_H1]={} MPa,[σ_H2]={} MPa'.format(simga_H_avilable_1,simga_H_avilable_2))
#     print('\t        [σ_H]=min([σ_H1],[σ_H2])={} MPa'.format(sigma_H_avilable))

#     print('\t   2)试算小齿轮分度圆直径')
#     d_1t_min=((2*K_Ht*T1/phi_d)*(1+1/u)*(Z_H*Z_E*Z_epsilon*Z_beta/sigma_H_avilable)**2)**(1/3)
#     print('\t     d_1t>={}mm'.format(d_1t_min))

#     print('\t(2)调整小齿轮分度圆直径')
#     print('\t   1)计算实际载荷系数前的数据准备')
#     v=pi*d_1t_min*n1/60/1000
#     b=phi_d*d_1t_min
#     print('\t     圆周速度v={}mm/s'.format(v))
#     print('\t     齿宽b={}mm'.format(v))
#     print('\t   2)计算实际载荷系数K_H')
#     K_A=1.25#使用系数，查表10-2取
#     K_V=1#查图10-8取
#     F_t1=2*T1/d_1t_min
#     if K_A*F_t1/b>=100:
#         K_H_alpha=1.2#对于7级精度的选择
#     else:
#         K_H_alpha=1.4#
#     K_H_beta=1.320#查表10-34得到
#     K_H=K_A*K_V*K_H_alpha*K_H_beta
#     d_1H=d_1t_min*(K_H/K_Ht)**(1/3)
#     m_nH=d_1H*cos(radians(beta))/z1
#     print('\t     使用系数 KA=',K_A)
#     print('\t     动载系数 KV=',K_V)
#     print('\t     计算齿间载荷分布系数 K_Hα=')
#     print('\t         小齿轮的径向力 F_t1=',F_t1)
#     print('\t         齿间载荷分布系数 K_Hα=',K_H_alpha)
#     print('\t     实际使用系数 K_H=',K_H)
#     print('\t     使用实际载荷系算得的分度圆直径 d_1H={}mm'.format(d_1H))
#     print('\t     相应的齿轮模数 m_nH,',m_nH)

#     print('\n3.按齿面弯曲疲劳强度设计')
#     print('\t(1)按式(10-17)试算模数')
#     print('\t   1)确定公式中的各种参数')
#     K_Ft=1.3
#     epsilon_alpha_v=epsilon_alpha/(cos(beta_b))**2
#     Y_epsilon=0.25+0.75/epsilon_alpha_v
#     Y_beta=1-epsilon_beta*beta/120
#     z_v1=z1/(cos(radians(beta)))**3
#     z_v2=z2/(cos(radians(beta)))**3
#     Y_Fa1=2.76#查表10-5得到
#     Y_Fa2=2.21#查表10-5得到
#     Y_Sa1=1.56#查表10-5得到
#     Y_Sa2=1.77#查表10-5得到
#     sigma_Flim1=1050#查图10-20c得
#     sigma_Flim2=1050#查图10-20c得
#     K_FN1=0.98#查图10-18得
#     K_FN2=1.05#查图10-18得
#     S=1.4#弯曲疲劳安全系数
#     sigma_F_avaliable_1=K_FN1*sigma_Flim1/S
#     sigma_F_avaliable_2=K_FN2*sigma_Flim2/S
#     num1=Y_Fa1*Y_Sa1/sigma_F_avaliable_1
#     num2=Y_Fa2*Y_Sa2/sigma_F_avaliable_2
#     num=max(num1,num2)
#     print('\t     试选载荷系数 K_Ft=',K_Ft)
#     print('\t     计算重合度系数 Y_ε') 
#     print('\t         ε_αv=',epsilon_alpha_v)
#     print('\t         重合度系数 Y_ε=',Y_epsilon) 
#     print('\t     螺旋角系数 Y_β=',Y_beta)
#     print('\t     计算(Y_Fa*Y_Sa/[σ_F])')
#     print('\t         当量齿数 Z_v1={},Z_v2={}'.format(z_v1,z_v2))
#     print('\t         齿形修正系数 Y_Fa1={},Y_Fa2={}'.format(Y_Fa1,Y_Fa2)) 
#     print('\t         应力修正系数 Y_Sa1={},Y_Sa2={}'.format(Y_Sa1,Y_Sa2)) 
#     print('\t         弯曲疲劳寿命系数 K_FN1={},K_FN2={}'.format(Y_Sa1,Y_Sa2))
#     print('\t         弯曲疲劳安全系数 S=',S)
#     print('\t         弯曲应力许用值 [σ_F]_1={},[σ_F]_2={}'.format(sigma_F_avaliable_1,sigma_F_avaliable_2))
#     print('\t         Y_Fa1*Y_Sa1/[σ_F1]={},Y_Fa2*Y_Sa2/[σ_F2]={}'.format(num1,num2))
#     print('\t         Y_Fa*Y_Sa/[σ_F]=max(Y_Fa1*Y_Sa1/[σ_F1],Y_Fa2*Y_Sa2/[σ_F2])={}'.format(num))
#     print('\t   2)试算齿轮模数')
#     m_nt_min=((2*K_Ft*T1*Y_epsilon*Y_beta*(cos(radians(beta)))**2)/(phi_d*z1**2)*num)**(1/3)
#     print('\t    m_nt_min=',m_nt_min)
#     print('\t(2)调整齿轮模数')
#     print('\t   1)计算实际载荷系数前的数据准备')
#     d1=m_nt_min*z1/cos(radians(beta))
#     v=pi*d1*n1/60/1000
#     b=phi_d*d1
#     h=(2*h_an_star+c_n_star)*m_nt_min
#     aspect_ratio=b/h
#     print('\t     计算圆周速度 v')
#     print('\t         小齿轮分度圆直径 d1=',d1)
#     print('\t         圆周速度 v=',v) 
#     print('\t     齿宽 b=',b)
#     print('\t     宽高比 b/h=',aspect_ratio)
#     print('\t   2)计算实际载荷系数')
#     K_A=1
#     K_V=1
#     F_t1=2*T1/d1
#     if K_A*F_t1/b>100:
#         K_F_alpha=1.2
#     else:
#         K_F_alpha=1.4
#     K_H_beta=1.250#插值法，查表10-4得
#     K_F_beta=1.20#查图 10-13得
#     K_F=K_A*K_V*K_F_alpha*K_F_beta
#     print('\t     动载系数 K_A',K_A)
#     print('\t     计算齿间载荷分布系数 K_Fα')
#     print('\t         F_t1=',F_t1)
#     print('\t         K_A*F_t1/b =',K_A*F_t1/b)
#     print('\t         齿间载荷分布系数 K_Fα=',K_F_alpha)
#     print('\t     齿间载荷分布系数 K_Hβ={},K_Fβ={}'.format(K_H_beta,K_F_beta))
#     print('\t     实际载荷系数 K_F=',K_F)
#     print('\t   3)按实际载荷系数算得齿轮模数和小齿轮分度圆直径')
#     m_nF=m_nt_min*(K_F/K_Ft)**(1/3)
#     d_1F=m_nF*z1/cos(radians(beta))
#     m_n=ceil(max(m_nF,m_nH))
#     m_n=9
#     d1=min(d_1F,d_1H)
#     z10=d_1H*cos(radians(beta))/m_n
#     z1=ceil(z1)
#     z2=floor(z1*u)
#     print('\t     齿轮模数 m_nF={}mm'.format(m_nF))
#     print('\t     小齿轮分度圆直径 d_1F=',d_1F)
#     print('\t     对比按弯曲疲劳强度和接触疲劳强度的计算结果')
#     print('\t         圆整齿轮模数 m_n={}mm'.format(m_n))
#     print('\t         分度圆直径 d1=',d1)
#     print('\t     则小齿轮齿数 z1={}，取为 z1={}'.format(z10,z1))
#     print('\t     则大齿轮齿数 z2={}，取为 z2={}'.format(z1*u,z2))

#     print('\n4.几何尺寸计算')
#     print('\t(1)计算中心距')
#     a0=(z1+z2)*m_n/2/cos(radians(beta))
#     a=ceil(a0)
#     print('\t      计算得中心距 a=',a0)
#     print('\t      圆整得中心距 a=',a)
#     print('\t(2)按圆整后的中心距修正螺旋角')
#     beta=degrees(acos((z1+z2)*m_n/2/a))
#     print('\t      修正后的螺旋角 β={}°'.format(beta))
#     print('\t(3)计算大小齿轮的分度圆直径')
#     d1=z1*m_n/cos(radians(beta))
#     d2=z2*m_n/cos(radians(beta))
#     print('\t      小齿轮分度圆直径 d1={}'.format(d1))
#     print('\t      大齿轮分度圆直径 d2={}'.format(d2))
#     print('\t(4)计算齿轮宽度')
#     b=phi_d*d1
#     b1=ceil(b+5)
#     b2=ceil(b)
#     print('\t      设计齿轮宽度 b={}mm'.format(b))
#     print('\t      取小齿轮宽度 b1={}mm'.format(b1))
#     print('\t      取大齿轮宽度 b2={}mm'.format(b2))
