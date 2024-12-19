
import numpy as np
from ngsolve.la import InnerProduct
from math import sqrt
from ngsolve import Variation
from ngsolve import grad
from ngsolve import TaskManager
from ngsolve import solvers
import scipy.sparse as sp
import time
# import matplotlib.pyplot as plt

COUNT=0

def NumpyNewton(u,dF,ddF,F=0,tol=1e-10,maxit=100,backtracking=True,damping=1):
    if F==0:
        backtracking=False
    for i in range(maxit):
        print("Newton iteration {:3}".format(i))
        print(" State "+str(u))
        print(" Energy "+str(F(u)))
        print(" Residual "+str(dF(u)))
        
        du=np.linalg.solve(ddF(u),-dF(u))
        m=np.dot(du,dF(u))
        if np.absolute(m)<tol:
            break
        alpha=1
        tau=0.5
        c=0.5
        t=-c*m
        if backtracking:
            for j in range(maxit):
                if F(u)-F(u+alpha*du)>=alpha*t:
                    break
                alpha=tau*alpha
                if j==maxit-1:
                    print("no descending stepsize found")
        else:
            alpha=min(1,damping*(i+1))
        u=u+alpha*du
        
        
def NGSLoadstepping(a,u,rhs,loadsteps=[1],hasenergy=False,tol=1e-16,maxit=100,backtracking=True,damping=1,inverse="pardiso",prnt=True):
    for l,ls in enumerate(loadsteps):
        if not loadsteps==[1]:
            input("Loadstep {:4}".format(l))   
        a+=Variation(ls*rhs)
        if not l==0:
            a+=Variation(-loadsteps[l-1]*rhs)
        NGSNewton(a,u,hasenergy,tol,maxit,backtracking,damping,inverse,prnt)

def NGSNewtonROTdc(a,B,Bt,D,Dt,u,l,hasenergy=False,tol=1e-16,maxit=100,backtracking=True,damping=1,inverse="pardiso",prnt=True,simplified=False):
    global COUNT
    if simplified:
        alphaBool=[False]*10
    else:
        alphaBool=[False]*2
    alphaCount=0
    if prnt: print("##### Newton start ######")
    Vn=B.shape[1]
    Xn=B.shape[0]+Vn
    Zn=D.shape[0]+Xn
    npu=np.zeros(Zn)
    npdu=np.zeros(Zn)
    npuh=np.zeros(Zn)
    npres=np.zeros(Zn)
    npresh=np.zeros(Zn)
    
    uh=u.vec.CreateVector()
    res=u.vec.CreateVector()
    resh=u.vec.CreateVector()
    npu[0:Vn]=np.array(u.vec)
    
        
    for i in range(maxit):
        if prnt: print("Newton iteration {:3}".format(i))
        # print("apply")
        a.Apply(u.vec,res)
        npres[0:Vn]=np.array(res)
        # print("assemlbe")
        if simplified: #Assemlbe only in First iteration
            if i<1:
                a.AssembleLinearization(u.vec)
        else:
            a.AssembleLinearization(u.vec)
        rows,cols,vals = a.mat.COO()  
        A = sp.csr_matrix((vals,(rows,cols)))   
        npres[0:Vn]+=Bt*npu[Vn:Xn]+Dt*npu[Xn:Zn]
        npres[Vn:Xn]=B*npu[0:Vn]
        npres[Xn:Zn]=D*npu[0:Vn]
        M=sp.bmat([[A,Bt,Dt],[B,None,None],[D,None,None]],format='csr')
        npdu=-sp.linalg.spsolve(M,npres)
        COUNT+=1
        
        # M=BlockMatrix( [[a.mat,bR.mat.T],[bR.mat,None]])
        # rhs=BlockVector( [-res-bR.mat.T*lam.vec,-bR.mat*u.vec])
        # preB=-bR.mat @ C.mat @ bR.mat.T
        # D=BlockMatrix([[C.mat,None],[None,preB]])
        # solvers.MinRes(mat=M,pre=D,rhs=rhs,sol=sol,tol=1e-3,maxsteps=1000,printrates=True)
        
        # du.data = -(a.mat+bR.mat).Inverse(freedofs=u.space.FreeDofs(),inverse=inverse)*res
        # print("solvedone")
        m=np.dot(npres,npdu)
        if hasenergy:
            energy=a.Energy(u.vec)+np.dot(npu[Vn:Xn],B*npu[0:Vn])+np.dot(npu[Xn:Zn],D*npu[0:Vn])
            if prnt: print(" Energy ", energy)
        else:
            # for fd in fixds:
            #     res[fd]=0
            energy=np.dot(npres,npres)
            if prnt: print(" Res ", energy)
        
        if prnt: print(" Error ", sqrt(abs(m)))
        if abs(m)<tol:
            print(i," Newtonsteps")
            # u.vec.data += du
            # gfuN.vec.data+=duN
            break
        tau=0.5
        alpha=1
        # print("ls")
        if backtracking:
            if hasenergy:
                c=0.499
                t=-c*m
                for j in range(min(maxit,20)):
                    npuh=npu+alpha*npdu
                    uh.FV().NumPy()[:]=npuh[0:Vn]
                    if energy-a.Energy(uh)-np.dot(npuh[Vn:Xn],B*npuh[0:Vn])-np.dot(npuh[Xn:Zn],D*npuh[0:Vn])>=alpha*t:
                        break
                    alpha=tau*alpha
                if prnt: print(" Linesearchsteps",j)
            else:
                mu=0.1
                for j in range(min(maxit,20)):
                    npuh=npu+alpha*npdu
                    uh.FV().NumPy()[:]=npuh[0:Vn]
                    a.Apply(uh,resh)
                    npresh[0:Vn]=np.array(resh)
                    # for fd in fixds:
                    #     resh[fd]=0    
                    npresh[0:Vn]+=Bt*npuh[Vn:Xn]+Dt*npuh[Xn:Zn]
                    npresh[Vn:Xn]=B*npuh[0:Vn]
                    npresh[Xn:Zn]=D*npuh[0:Vn]
                    if np.dot(npresh,npresh)<=energy:
                        break
                    alpha=tau*alpha
                if prnt: print(" Linesearchsteps",j)
                
        else:
            alpha=min(damping*(i+1),1)
        if backtracking:
            if alpha==1:
                alphaBool[alphaCount]=True
                alphaCount+=1
            if np.prod(alphaBool):
                backtracking=False
            if (j==min(maxit,20)-1):
                print("no descending stepsize found")
                break   
        npu += alpha*npdu
        u.vec.FV().NumPy()[:]=npu[0:Vn]
        l[:]=npu[Vn:Xn]


def NGSNewtonROT2(a,B,Bt,u,l,hasenergy=False,tol=1e-16,maxit=100,backtracking=True,damping=1,inverse="pardiso",prnt=True,simplified=False):
    if simplified:
        alphaBool=[False]*10
    else:
        alphaBool=[False]*2
    alphaCount=0
    if prnt: print("##### Newton start ######")
    Vn=B.shape[1]
    Xn=B.shape[0]+Vn
    npu=np.zeros(Xn)
    npdu=np.zeros(Xn)
    npuh=np.zeros(Xn)
    npres=np.zeros(Xn)
    npresh=np.zeros(Xn)
    
    uh=u.vec.CreateVector()
    res=u.vec.CreateVector()
    resh=u.vec.CreateVector()
    npu[0:Vn]=np.array(u.vec)
    
        
    for i in range(maxit):
        if prnt: print("Newton iteration {:3}".format(i))
        # print("apply")
        a.Apply(u.vec,res)
        npres[0:Vn]=np.array(res)
        # print("assemlbe")
        if simplified: #Assemlbe only in First iteration
            if i<1:
                a.AssembleLinearization(u.vec)
        else:
            a.AssembleLinearization(u.vec)
        rows,cols,vals = a.mat.COO()  
        A = sp.csr_matrix((vals,(rows,cols)))   
        npres[0:Vn]+=Bt*npu[Vn:Xn]
        npres[Vn:Xn]=B*npu[0:Vn]
        M=sp.bmat([[A,Bt],[B,None]],format='csr')
        npdu=-sp.linalg.spsolve(M,npres)
        
        # M=BlockMatrix( [[a.mat,bR.mat.T],[bR.mat,None]])
        # rhs=BlockVector( [-res-bR.mat.T*lam.vec,-bR.mat*u.vec])
        # preB=-bR.mat @ C.mat @ bR.mat.T
        # D=BlockMatrix([[C.mat,None],[None,preB]])
        # solvers.MinRes(mat=M,pre=D,rhs=rhs,sol=sol,tol=1e-3,maxsteps=1000,printrates=True)
        
        # du.data = -(a.mat+bR.mat).Inverse(freedofs=u.space.FreeDofs(),inverse=inverse)*res
        # print("solvedone")
        m=np.dot(npres,npdu)
        if hasenergy:
            energy=a.Energy(u.vec)+np.dot(npu[Vn:Xn],B*npu[0:Vn])
            if prnt: print(" Energy ", energy)
        else:
            # for fd in fixds:
            #     res[fd]=0
            energy=np.dot(npres,npres)
            if prnt: print(" Res ", energy)
        
        if prnt: print(" Error ", sqrt(abs(m)))
        if abs(m)<tol:
            print(i," Newtonsteps")
            # u.vec.data += du
            # gfuN.vec.data+=duN
            break
        tau=0.5
        alpha=1
        # print("ls")
        if backtracking:
            if hasenergy:
                c=0.499
                t=-c*m
                for j in range(min(maxit,20)):
                    npuh=npu+alpha*npdu
                    uh.FV().NumPy()[:]=npuh[0:Vn]
                    if energy-a.Energy(uh)-np.dot(npuh[Vn:Xn],B*npuh[0:Vn])>=alpha*t:
                        break
                    alpha=tau*alpha
                if prnt: print(" Linesearchsteps",j)
            else:
                mu=0.1
                for j in range(min(maxit,20)):
                    npuh=npu+alpha*npdu
                    uh.FV().NumPy()[:]=npuh[0:Vn]
                    a.Apply(uh,resh)
                    npresh[0:Vn]=np.array(resh)
                    # for fd in fixds:
                    #     resh[fd]=0    
                    npresh[0:Vn]+=Bt*npuh[Vn:Xn]
                    npresh[Vn:Xn]=B*npuh[0:Vn]
                    if np.dot(npresh,npresh)<=energy:
                        break
                    alpha=tau*alpha
                if prnt: print(" Linesearchsteps",j)
                
        else:
            alpha=min(damping*(i+1),1)
        if backtracking:
            if alpha==1:
                alphaBool[alphaCount]=True
                alphaCount+=1
            if np.prod(alphaBool):
                backtracking=False
            if (j==min(maxit,20)-1):
                print("no descending stepsize found")
                break   
        npu += alpha*npdu
        u.vec.FV().NumPy()[:]=npu[0:Vn]
        l[:]=npu[Vn:Xn]
        
        
def NGSNewtonROT(a,bR,u,hasenergy=False,tol=1e-16,maxit=100,backtracking=True,damping=1,inverse="pardiso",prnt=True,simplified=False):
    if simplified:
        alphaBool=[False]*10
    else:
        alphaBool=[False]*2
    alphaCount=0
    fixds=list(np.where(np.array(u.space.FreeDofs())*1<0.1)[0])
    if prnt: print("##### Newton start ######")
    du=u.vec.CreateVector()
    res=u.vec.CreateVector()
    uh=u.vec.CreateVector()
    resh=u.vec.CreateVector()
    Brows,Bcols,Bvals = bR.mat.COO()  
    B = sp.csr_matrix((Bvals,(Brows,Bcols)))
    C = sp.csr_matrix((Bvals,(Bcols,Brows)))
    Vn=B.shape[1]
    Xn=B.shape[0]+Vn
    
    for i in range(maxit):
        if prnt: print("Newton iteration {:3}".format(i))
        # print("apply")
        a.Apply(u.vec[0:Vn],res[0:Vn])
        # print("assemlbe")
        if simplified: #Assemlbe only in First iteration
            if i<1:
                a.AssembleLinearization(u.vec[0:Vn])
        else:
            a.AssembleLinearization(u.vec[0:Vn])
        rows,cols,vals = a.mat.COO()  
        A = sp.csr_matrix((vals,(rows,cols)))          
        res.data[0:Vn]+=bR.mat.T*u.vec[Vn:Xn]
        res.data[Vn:Xn]=bR.mat*u.vec[0:Vn]
        M=sp.bmat([[A,C],[B,None]],format='csr')
        du.FV().NumPy()[:]=-sp.linalg.spsolve(M,res.FV().NumPy()[:])
        
        # M=BlockMatrix( [[a.mat,bR.mat.T],[bR.mat,None]])
        # rhs=BlockVector( [-res-bR.mat.T*lam.vec,-bR.mat*u.vec])
        # preB=-bR.mat @ C.mat @ bR.mat.T
        # D=BlockMatrix([[C.mat,None],[None,preB]])
        # solvers.MinRes(mat=M,pre=D,rhs=rhs,sol=sol,tol=1e-3,maxsteps=1000,printrates=True)
        
        # du.data = -(a.mat+bR.mat).Inverse(freedofs=u.space.FreeDofs(),inverse=inverse)*res
        # print("solvedone")
        m=InnerProduct(res,du)
        if hasenergy:
            energy=a.Energy(u.vec[0:Vn])+InnerProduct(u.vec[Vn:Xn],bR.mat*u.vec[0:Vn])
            if prnt: print(" Energy ", energy)
        else:
            # for fd in fixds:
            #     res[fd]=0
            energy=InnerProduct(res,res)
            if prnt: print(" Res ", energy)
        
        if prnt: print(" Error ", sqrt(abs(m)))
        if abs(m)<tol:
            print(i," Newtonsteps")
            # u.vec.data += du
            # gfuN.vec.data+=duN
            break
        tau=0.5
        alpha=1
        # print("ls")
        if backtracking:
            if hasenergy:
                c=0.499
                t=-c*m
                for j in range(min(maxit,20)):
                    uh.data=u.vec+alpha*du
                    if energy-a.Energy(uh[0:Vn])-InnerProduct(uh[Vn:Xn],bR.mat*uh[0:Vn])>=alpha*t:
                        break
                    alpha=tau*alpha
                if prnt: print(" Linesearchsteps",j)
            else:
                mu=0.1
                for j in range(min(maxit,20)):
                    uh.data=u.vec+alpha*du
                    a.Apply(uh[0:Vn],resh[0:Vn])
                    # for fd in fixds:
                    #     resh[fd]=0    
                    resh.data[0:Vn]+=bR.mat.T*uh[Vn:Xn]
                    resh.data[Vn:Xn]=bR.mat*uh[0:Vn]
                    if InnerProduct(resh,resh)<=energy:
                        break
                    alpha=tau*alpha
                if prnt: print(" Linesearchsteps",j)
                
        else:
            alpha=min(damping*(i+1),1)
        if backtracking:
            if alpha==1:
                alphaBool[alphaCount]=True
                alphaCount+=1
            if np.prod(alphaBool):
                backtracking=False
            if (j==min(maxit,20)-1):
                print("no descending stepsize found")
                break   
        u.vec.data += alpha*du  
        
        

def NGSNewton(a,u,C=None,hasenergy=False,tol=1e-16,maxit=100,backtracking=True,damping=1,inverse="pardiso",prnt=True,simplified=False):
    if simplified:
        alphaBool=[False]*10
    else:
        alphaBool=[False]*2
    alphaCount=0
    fixds=list(np.where(np.array(u.space.FreeDofs())*1<0.1)[0])
    if prnt: print("##### Newton start ######")
    du=u.vec.CreateVector()
    res=u.vec.CreateVector()
    uh=u.vec.CreateVector()
    resh=u.vec.CreateVector()
    
    for i in range(maxit):
        with TaskManager():
            # t1=time.time()
            if prnt: print("Newton iteration {:3}".format(i))
            a.Apply(u.vec,res)
            # t2=time.time()
            # print("apply",t2-t1)
            if simplified:
                if i<1:
                    a.AssembleLinearization(u.vec)
            else:
                a.AssembleLinearization(u.vec)
            # t3=time.time()
            # print("assemble",t3-t2)
            if C==None:
                du.data = -a.mat.Inverse(freedofs=u.space.FreeDofs(),inverse=inverse)*res
            else:
                solvers.MinRes(mat=a.mat,pre=C,rhs=-res,sol=du,printrates=True)
            # t4=time.time()
            # print("solve",t4-t3)
            m=InnerProduct(du,res)
            if hasenergy:
                energy=a.Energy(u.vec)
                if prnt: print(" Energy ", energy)
            else:
                for fd in fixds:
                    res[fd]=0
                energy=InnerProduct(res,res)
                if prnt: print(" Res ", energy)
            
            if prnt: print(" Error ", sqrt(abs(m)))
            if abs(m)<tol:
                print(i," Newtonsteps")
                break
            tau=0.5
            alpha=1
            if backtracking:
                if hasenergy:
                    c=0.499
                    t=-c*m
                    for j in range(min(maxit,20)):
                        uh.data=u.vec+alpha*du
                        if energy-a.Energy(uh)>=alpha*t:
                            break
                        alpha=tau*alpha
                    if prnt: print(" Linesearchsteps",j)
                else:
                    mu=0.1
                    for j in range(min(maxit,20)):
                        uh.data=u.vec+alpha*du
                        a.Apply(uh,resh)
                        for fd in fixds:
                            resh[fd]=0                    
                        if InnerProduct(resh,resh)<=energy:
                            break
                        alpha=tau*alpha
                    if prnt: print(" Linesearchsteps",j)
                    
            else:
                alpha=min(damping*(i+1),1)
            if backtracking:
                if alpha==1:
                    alphaBool[alphaCount]=True
                    alphaCount+=1
                if np.prod(alphaBool):
                    backtracking=False
                if (j==min(maxit,20)-1):
                    print("no descending stepsize found")
                    break    
            # t5=time.time()
            # print("linesearch",t5-t4)
        u.vec.data += alpha*du

def NGSNewtonPre(a,u,C,hasenergy=False,tol=1e-16,maxit=100,backtracking=True,damping=1,inverse="pardiso",prnt=True,simplified=False):
    if simplified:
        alphaBool=[False]*10
    else:
        alphaBool=[False]*2
    alphaCount=0
    fixds=list(np.where(np.array(u.space.FreeDofs())*1<0.1)[0])
    if prnt: print("##### Newton start ######")
    print("ndof",u.space.ndof)
    du=u.vec.CreateVector()
    res=u.vec.CreateVector()
    uh=u.vec.CreateVector()
    resh=u.vec.CreateVector()
    
    for i in range(maxit):
        with TaskManager():
            t1=time.time()
            if prnt: print("Newton iteration {:3}".format(i))
            a.Apply(u.vec,res)
            t2=time.time()
            print("apply",t2-t1)
            if simplified:
                if i<1:
                    a.AssembleLinearization(u.vec)
            else:
                a.AssembleLinearization(u.vec)
            t3=time.time()
            print("assemble",t3-t2)
            # if i<1:
            #     C=a.mat.Inverse(u.space.FreeDofs(),inverse=inverse)    
            #     print("preconditioner",time.time()-t3)
            #     t3=time.time()
            solvers.GMRes(a.mat,-res,pre=C.mat,x=du,printrates=True,maxsteps=1000)
            # du=solvers.PreconditionedRichardson(a,-res,pre=C,printing=True,maxit=1000)
            t4=time.time()
            print("solve",t4-t3)
            m=InnerProduct(du,res)
            if hasenergy:
                energy=a.Energy(u.vec)
                if prnt: print(" Energy ", energy)
            else:
                for fd in fixds:
                    res[fd]=0
                energy=InnerProduct(res,res)
                if prnt: print(" Res ", energy)
            
            if prnt: print(" Error ", sqrt(abs(m)))
            if abs(m)<tol:
                print(i," Newtonsteps")
                break
            tau=0.5
            alpha=1
            if backtracking:
                if hasenergy:
                    c=0.499
                    t=-c*m
                    for j in range(min(maxit,20)):
                        uh.data=u.vec+alpha*du
                        if energy-a.Energy(uh)>=alpha*t:
                            break
                        alpha=tau*alpha
                    if prnt: print(" Linesearchsteps",j)
                else:
                    mu=0.1
                    for j in range(min(maxit,20)):
                        uh.data=u.vec+alpha*du
                        a.Apply(uh,resh)
                        for fd in fixds:
                            resh[fd]=0                    
                        if InnerProduct(resh,resh)<=energy:
                            break
                        alpha=tau*alpha
                    if prnt: print(" Linesearchsteps",j)
                    
            else:
                alpha=min(damping*(i+1),1)
            if backtracking:
                if alpha==1:
                    alphaBool[alphaCount]=True
                    alphaCount+=1
                if np.prod(alphaBool):
                    backtracking=False
                if (j==min(maxit,20)-1):
                    print("no descending stepsize found")
                    break    
            t5=time.time()
            print("linesearch",t5-t4)
        u.vec.data += alpha*du

def NGSGFNewton(a,u,f,rhs,F,B,RotM,loadsteps=[1],damping=False,tol=1e-16,maxit=100,inverse="pardiso"):
    for l,ls in enumerate(loadsteps):
        if not loadsteps==[1]:
            print("Loadstep {:4}".format(l))    
        f+=ls*rhs
        if not l==0:
            f+=-loadsteps[l-1]*rhs
        res=u.vec.CreateVector()
        du=u.vec.CreateVector()
        uh=u.vec.CreateVector()
        resh=u.vec.CreateVector()
        for i in range(maxit):
            print("Newton iteration {:3}".format(i))          
            B.Set(RotM*grad(u))
            F(np.array(B.vec))
            a.Assemble()
            f.Assemble()
            du.data = a.mat.Inverse(u.space.FreeDofs(),inverse=inverse)*f.vec
            res.data=-f.vec
            m=InnerProduct(du,res)
            print(" Error ", sqrt(abs(m)))
            if abs(m)<tol:
                break
            alpha=1
            if damping:
                tau=0.5
                mu=0.1
                for j in range(maxit):
                    if j==0:
                        u.vec.data += alpha*du
                    else:
                        u.vec.data += (1-1/tau)*alpha*du
                    B.Set(RotM*grad(u))
                    F(np.array(B.vec))
                    f.Assemble
                    if InnerProduct(du,f.vec)<=-(1-2*mu*alpha)*m:
                        break
                    alpha=tau*alpha
                    if j==maxit-1:
                        input("no descending stepsize found")
                print(" Dampingsteps",j)
            else:
                u.vec.data += alpha*du