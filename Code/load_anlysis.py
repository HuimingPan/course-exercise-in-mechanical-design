from math import tan,cos,radians

def force_analyze(gear):
    alpha=radians(gear.alpha)
    beta=radians(gear.beta)
    F_t1=2*gear.T/gear.d1*1000
    F_r1=F_t1*tan(alpha)/cos(beta)
    F_a1=F_t1*tan(beta)
    print("F_t1={:.7} N,F_r1={:.7} N,F_a1={:.7} N".format(F_t1,F_r1,F_a1))