# 此脚本用于找Rys多项式的递推系数

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def integ(f1,f2)->float:
    """计算内积 ∫f1(x)f2(x)w(x)dx"""
    func=lambda x: f1(x)*f2(x)*wfunc(x)
    return quad(func,up,lo)[0]

def alp_func(f)->float:
    """计算 α_k = 〈x p_k, p_k〉 / 〈p_k, p_k〉"""
    xp=lambda x:x*f(x)
    return integ(xp,f)/integ(f,f)

def bet_func(f,fn)->float:
    """计算 β_k = 〈p_k, p_k〉 / 〈p_{k-1}, p_{k-1}〉"""
    return integ(f,f)/integ(fn,fn)


alphas=np.linspace(0,1,50)
alphas=alphas**(1.5)*30

roots=[]
alpl=[]
betl=[]
N=20
up,lo=-1.0,1.0
# 使用平台支持的最长浮点类型
float_type = np.longdouble
print(f"使用浮点类型: {float_type}")
for alpha in alphas:

    wfunc = lambda x: np.exp(-alpha*x**2) # 权函数

    _alps=np.zeros(N,dtype=float_type) # 递推系数
    _bets=np.zeros(N,dtype=float_type)


    alps=lambda i:_alps[i]
    bets=lambda i:_bets[i]

    p0=lambda x: 1.0
    p1=lambda x: x-alps(0)

    ps=[p0,p1]
    namespace = {'alps': alps, 'bets': bets, 'p0': p0, 'p1': p1}
    for i in range(2,N+1):
        form=f"lambda x: (x-alps({i-1}))*p{i-1}(x)-bets({i-1})*p{i-2}(x)"
        # print(i,form)
        p_new = eval(form, namespace)
        ps.append(p_new)
        # 将新创建的函数添加到命名空间中供后续函数使用
        namespace[f'p{i}'] = p_new
    
    _bets[0]=integ(p0,p0)
    for i in range(N):
        _alps[i]=alp_func(ps[i])
        if i==0:continue
        _bets[i]=bet_func(ps[i],ps[i-1])
    alpl.append(_alps)
    betl.append(_bets)

alpl=np.array(alpl)
betl=np.array(betl)
print(alpl.shape)
print(betl.shape)

np.savetxt("betl.txt",betl.T) # 存储递推系数 β，α都是0(为啥捏)
np.savetxt("alps.txt",alphas)