import pyUni10 as uni10
import sys
import numpy as np

def matId():
  spin = 1.0
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [1, 0, 0, 0, 1, 0, 0, 0, 1])

def matSm():
  spin = 1.0
  dim = int(spin * 2 + 1)
  return np.sqrt(2)*uni10.Matrix(dim, dim, [0, 0, 0, 1, 0, 0, 0, 1, 0])

def matSp():
  spin = 1.0
  dim = int(spin * 2 + 1)
  return np.sqrt(2)*uni10.Matrix(dim, dim, [0, 1, 0, 0, 0, 1, 0, 0, 0])

def matSz():
  spin = 1.0
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [1, 0, 0, 0, 0, 0, 0, 0, -1])

def matSz2():
  spin = 1.0
  dim = int(spin * 2 + 1)
  return uni10.Matrix(dim, dim, [1, 0, 0, 0, 0, 0, 0, 0, 1])

def Heisenberg(D, Delta, h, J):
    iden = matId()
    sm = matSm()
    sp = matSp()
    sz = matSz()
    sz2 = matSz2()
    ham = float(J)*(0.5*(uni10.otimes(sm,sp)+uni10.otimes(sp,sm))+float(Delta)*uni10.otimes(sz,sz))
    ham += 0.5*float(h)*(uni10.otimes(iden,sz)+uni10.otimes(sz,iden))
    ham += 0.5*float(D)*(uni10.otimes(iden,sz2)+uni10.otimes(sz2,iden))
    spin = 1.0
    dim = int(spin * 2 + 1)
    bdi = uni10.Bond(uni10.BD_IN, dim);
    bdo = uni10.Bond(uni10.BD_OUT, dim);
    H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
    H.putBlock(ham)
    return H

def bondcat(T, L, bidx):
    labels = T.label();
    per_labels = list(T.label())
    per_labels.insert(0, per_labels.pop(bidx))
    inBondNum = T.inBondNum()
    T.permute(per_labels, 1)
    T.putBlock(L * T.getBlock())
    T.permute(labels, inBondNum)

def bondrm(T, L, bidx):
    invL = uni10.Matrix(L.row(), L.col(), True)
    for i in xrange(L.elemNum()):
        invL[i] = 0 if L[i] == 0 else 1.0 / L[i]
    bondcat(T, invL, bidx)
    
def findinv(chi, G):
    A = 1 * G
    B = 1 * G
    B.transpose()
    A.setLabel([1, 3, 0])
    B.setLabel([2, 0, 4])
    T = A * B
    T.permute([1, 2, 3, 4], 2)
    #E=T.getBlock()
    #F=T.getBlock().transpose()
    #for i in range(chi*chi*chi*chi):
    #    if (E[i]!=F[i]):
    #        print E[i],F[i]
    D, U = T.getBlock().eigh()
    lambda1 = 0.000
    for i in range(chi*chi):
        if (abs(D[i])>lambda1):
            lambda1 = abs(D[i])
            i1 = i
    return D[i1]
    
def findxin(chi, Ntheta):
    Ntheta.setLabel([1, 2, 3, 4])
    Mtheta = 1 * Ntheta
    Mtheta.setLabel([-1, 2, -3, 4])
    E = Mtheta * Ntheta
    E.permute([1, -1, 3, -3], 2)
    D, U = E.getBlock().eigh()
    lambda1 = 0.000
    lambda2 = 0.000
    for i in range(chi*chi):
        if (abs(D[i])>lambda1):
            lambda1 = abs(D[i])
            i1 = i
    for i in range(chi*chi):
        if i!=i1:
            if (abs(D[i])>lambda2):
                lambda2 = abs(D[i])
    return -1/np.log(np.sqrt(lambda2/lambda1))

def findxis(chi, A, B):
    A.setLabel([1, 3, 0])
    B.setLabel([2, 4, 0])
    T = A * B
    T.permute([1, 2, 3, 4], 2)
    D, U = T.getBlock().eigh()
    lambda1 = 0.000
    lambda2 = 0.000
    i1=0
    for i in range(chi*chi):
        if (abs(D[i])>lambda1):
            lambda1 = abs(D[i])
            i1 = i
    for i in range(chi*chi):
        if i!=i1:
            if (abs(D[i])>lambda2):
                lambda2 = abs(D[i])
    return -1/np.log(lambda2/lambda1)

D = float(sys.argv[1])
chi = int(sys.argv[2])
istate = int(sys.argv[3])

delta = 0.01
N = 20000
Delta = 1.000
h = 0.000
J = 1.000
H = Heisenberg(D, Delta, h, J)

bdi_chi = uni10.Bond(uni10.BD_IN, chi)
bdo_chi = uni10.Bond(uni10.BD_OUT, chi)

# Gamma matrices: Gs=[Ga, Gb]
Gs = []
Gs.append(uni10.UniTensor([bdi_chi, bdo_chi, H.bond(2)], "Ga"))
Gs.append(uni10.UniTensor([bdi_chi, bdo_chi, H.bond(2)], "Gb"))
Gs[0].set_zero(), Gs[1].set_zero()
#Gs[0].randomize(), Gs[1].randomize()
#Gs.append(uni10.UniTensor([bdi_chi, bdo_chi, H.bond(2)], "Gc"))
#Gs.append(uni10.UniTensor([bdi_chi, bdo_chi, H.bond(2)], "Gd"))
#Gs[2].set_zero(), Gs[3].set_zero()

# Lambda matrices (diagonal): Ls=[La, Lb]
Ls = []
Ls.append(uni10.Matrix(chi, chi, True))  # Diagonal matrix
Ls.append(uni10.Matrix(chi, chi, True))  # Diagonal matrix
Ls[0].set_zero(), Ls[1].set_zero()
#Ls[0].randomize(), Ls[1].randomize()

# Setup U=exp^{-delta H}
U = uni10.UniTensor(H.bond(), "U")
U.putBlock(uni10.takeExp(-delta, H.getBlock()))

# Initialize a product state
if (istate != 2):
  Ls[0][0]=1.000
  Ls[1][0]=1.000 # D=1 Product State
else:
  Ls[0][0]=np.sqrt(1/2.)
  Ls[0][1]=np.sqrt(1/2.)
  Ls[1][0]=np.sqrt(1/2.)
  Ls[1][1]=np.sqrt(1/2.) # D=2 Product State

# Initialize a product state
Ms0 = Gs[0].getBlock()
Ms1 = Gs[1].getBlock()
if (istate == 0): # Large-D
  Ms0[1]=1.000 # Sz=0
  Ms1[1]=1.000 # Sz=0
elif (istate == 1): # Neel
  Ms0[0]=1.000 # Sz=+1
  Ms1[2]=1.000 # Sz=-1
elif (istate == 2): # AKLT
  Ms0[1]=-np.sqrt(2/3.) # A+
  Ms0[3]=np.sqrt(4/3.) # A0
  Ms0[2+3*chi]=-np.sqrt(4/3.) #A0
  Ms0[4+3*chi]=np.sqrt(2/3.) # A-
  Ms1[1]=-np.sqrt(2/3.) # A+
  Ms1[3]=np.sqrt(4/3.) # A0
  Ms1[2+3*chi]=-np.sqrt(4/3.) #A0
  Ms1[4+3*chi]=np.sqrt(2/3.) # A-
Gs[0].putBlock(Ms0)
Gs[1].putBlock(Ms1)

#Product state for strange correlator
#Gs[2].putBlock(Ms)

Elist=[]
#Ilist=[]
Mxlist=[]
Mzlist=[]
Slist=[]
#Xnlist=[]
#Xslist=[]
Ylist=[]

for level in range(chi):
    Ylist.append([])

for step in range(N):
    
    # Construct theta
    A = step % 2
    B = (step + 1) % 2
    bondcat(Gs[A], Ls[A], 1)
    bondcat(Gs[A], Ls[B], 0)
    bondcat(Gs[B], Ls[B], 1)
    #Gs[3] = 1 * Gs [B]
    Gs[A].setLabel([-1, 3, 1])
    Gs[B].setLabel([3, -3, 2])
    #theta = uni10.contract(Gs[A], Gs[B], False) # Gs[A], Gs[B] is NOT permuted after the execution
    theta = uni10.contract(Gs[A], Gs[B], True) # Gs[A], Gs[B] is permuted after the execution

    #bondrm(Gs[A], Ls[B], 0) # one or the other
    #bondrm(Gs[B], Ls[B], 1) # one or the other
    #Gs[A].setLabel([-1, 3, 1])
    #Gs[B].setLabel([3, -3, 2])
    #Ntheta = uni10.contract(Gs[A], Gs[B], True) # Gs[A], Gs[B] is permuted after the execution
    
    # Define magnetization
    Sx = 0.5*(matSm() + matSp()) #Mx
    Sz = matSz() #Mz
    Enorm = (theta * theta)[0]
    Mtheta = 1 * theta
    bondcat(Mtheta, Sx, 1)
    Mx = (Mtheta * theta)[0] / Enorm
    Mtheta = 1 * theta
    bondcat(Mtheta, Sz, 1)
    Mz = (Mtheta * theta)[0] / Enorm
    #Ntheta = 1 * theta

    # Evolve theta
    U.setLabel([1, 2, -2, -4])
    theta *= U;
    theta.permute([-1, -2, -3, -4], 2)

    # SVD
    svd = theta.getBlock().svd()

    # Truncation
    sv = svd[1]
    norm = sv.resize(chi, chi).norm()
    sv *= 1.0 / norm;
    Ls[A] = sv
    Gs[A].putBlock(svd[0].resize(svd[0].row(), chi))
    Gs[B].putBlock(svd[2].resize(chi, svd[2].col()))
    Gs[A].permute([-1, 3, 1], 1)
    bondrm(Gs[A], Ls[B], 0)
    bondrm(Gs[B], Ls[B], 1)  
    val = (theta * theta)[0]
    #norm =(Ntheta * Ntheta )[0]
    E = -np.log(val) / delta / 2
    E = E / Enorm
    #norm =(theta * theta )[0]
    #M = (Mtheta*theta)[0]/Enorm
    #M = (Mtheta * Ntheta)[0]/norm
    
    # Record Observables
    if step % 2 == 0:
    
        # Measure Entanglement
        S=0.00
        for level in range(chi):
            plevel=np.square(Ls[0][level])
            if (plevel>0):
                S=S-plevel*np.log(plevel)
                Ylist[level].append(-np.log(plevel))
            else:
                Ylist[level].append(0)
        
        #I=findinv(chi, Gs[3])
        #Xn=findxin(chi, Ntheta)
        #Xs=findxis(chi, Gs[2], Gs[3])
        
        Elist.append(E)
        #Ilist.append(I)
        Mxlist.append(Mx)
        Mzlist.append(Mz)
        Slist.append(S)
        #Xnlist.append(Xn)
        #Xslist.append(Xs)

f = open('entropy_D='+str(D)+'_chi='+str(chi)+'_istate='+str(istate)+'.dat', 'w')
for item in Slist:
  print >> f, item
f = open('energy_D='+str(D)+'_chi='+str(chi)+'_istate='+str(istate)+'.dat', 'w')
for item in Elist:
  print >> f, item
f = open('spectrum_D='+str(D)+'_chi='+str(chi)+'_istate='+str(istate)+'.dat', 'w')
for level in range(chi):
  for item in Ylist[level]:
    print >> f, item
  print >> f
