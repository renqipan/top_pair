# refer to: https://arxiv.org/abs/1305.1878
import numpy as np
import ROOT as r
import math
from scipy.optimize import leastsq

mT =172.5 #GeV : top quark mass
mW=80.385 #GeV: W boson mass
mN=0.0 #GeV: neutrino mass
def UnitCircle():
    '''Unit circle in extended representation'''
    return np.diag([1,1,-1])
def cofactor(A, ij):
    '''Cofactor[i,j] of 3x3 matrix A'''
    i, j=ij
    a=A[not i:2 if i==2 else None:2 if i==1 else 1,
        not j:2 if j==2 else None: 2 if j==1 else 1]
    return (-1)**(i+j)*(a[0,0]*a[1,1]-a[1,0]*a[0,1])
def cofactor2(A, ij):
    '''Cofactor[i,j] of 3x3 matrix A'''
    i, j=ij
    a = np.delete(np.delete(A, i, axis=0), j, axis=1)
    return (-1)**(i+j)*(a[0,0]*a[1,1]-a[1,0]*a[0,1])
def R(axis, angle):
    '''Rotation matrix about x(0),y(1),z(2) axis'''
    c,s=math.cos(angle),math.sin(angle)
    R=c*np.eye(3)
    for i in [-1,0,1]:
        R[(axis-i)%3,(axis+i)%3]=i*s+(1-i*i)
    return R
def Derivative():
    '''Matrix to differentiate[cos(t), sin(t),1]'''
    return  R(2,math.pi/2).dot(np.diag([1,1,0]))
def multisqrt(y):
    '''valide real solutions y=x*x'''
    return ([] if y< 0 else [0] if y==0 else 
            (lambda r:[-r, r])(math.sqrt(y)))
def factor_degenerate(G, zero=0):
    '''Linear factors of degenerate quadratic polynomial'''
    if G[0,0]==0==G[1,1]:
        return [[G[0,1],0,G[1,2]],[0,G[0,1],G[0,2]-G[1,2]]]
    swapXY=1 if abs(G[0,0])> abs(G[1,1]) else 0
    #swapXY=abs(G[0,0])> abs(G[1,1])
    Q=G[(1,0,2),][:,(1,0,2)] if swapXY else G
    Q=Q/Q[1,1]
    q22=cofactor(Q,(2,2))
    if -q22 <= zero:
        lines=[[Q[0,1],Q[1,1],Q[1,2]+s]
               for s in multisqrt(-cofactor(Q,(0,0)))]
    else:
        x0, y0=[cofactor(Q,(i,2))/q22 for i in [0,1]]
        lines=[[m,Q[1,1],-Q[1,1]*y0-m*x0] 
               for m in [Q[0,1]+s for s in multisqrt(-q22)]]
    return [[L[swapXY],L[not swapXY], L[2]] for L in lines]
def intersections_ellipse_line(ellipse,line, zero=1e-12):
    '''points of intersection between ellipse and line'''
    _,V=np.linalg.eig(np.cross(line,ellipse).T)
    sols=sorted([(v.real / v[2].real,
                  np.dot(line,v.real)**2+
                  np.dot(v.real,ellipse).dot(v.real)**2)
                for v in V.T], 
                 key=lambda aa:aa[1])[:2]
    return [s for s ,k in sols if k < zero]
def intersections_ellipses(A, B, returnLines=False):
    '''Points of intersection between two ellipse'''
    LA=np.linalg
    if abs(LA.det(B)) > abs(LA.det(A)): A, B =B, A
    e=next(e.real for e in LA.eigvals(LA.inv(A).dot(B))
           if not e.imag)
    lines=factor_degenerate(B- e*A)
    points=sum([intersections_ellipse_line(A,L) for L in lines], [])
    return (points, lines) if returnLines else points
class nuSolutionSet(object):
    '''Definition for nu analytic solution, t-> b, mu, nu'''
    def __init__(self, b, mu, #Lorentz vectors
                mW2=mW**2, mT2=mT**2, mN2=mN**2):
        c=r.Math.VectorUtil.CosTheta(b,mu)
        s=math.sqrt(1-c**2)

        x0p=-(mT2-mW2-b.M2())/(2*b.E())
        x0=-(mW2-mu.M2()-mN2)/(2*mu.E())
        Bb, Bm=b.Beta(), mu.Beta()
        Sx=(x0*Bm-mu.P()*(1-Bm**2))/Bm**2
        Sy=(x0p/Bb -c*Sx)/s
        w=(Bm / Bb-c) / s
        w_=(-Bm / Bb-c)/s
        Om2=w**2+1-Bm**2
        eps2=(mW2-mN2)*(1-Bm**2)
        x1=Sx - (Sx+w*Sy) / Om2
        y1=Sy-(Sx+w*Sy)*w / Om2
        Z2=x1**2*Om2-(Sy-w*Sx)**2-(mW2-x0**2-eps2)
        Z=math.sqrt(max(0, Z2))
        for item in ['b', 'mu','c','s','x0','x0p','Sx','Sy',
                    'w','w_','x1','y1','Z','Om2','eps2','mW2']:
            setattr(self, item, eval(item))
    @property
    def K(self):
        '''Extended rotation from F' to F coordinate'''
        return np.array([self.c,-self.s,0,0],[self.s,self.c,0,0],
                        [0,0,1,0],[0,0,0,1])
    @property
    def A_mu(self):
        '''F coordinate constraint on w momentum: ellipsoid'''
        B2=self.mu.Beta()**2
        SxB2=self.Sx*B2
        F=self.mW2-self.x0**2-self.eps2
        return np.array([[1-B2,0,0,SxB2],[0,1,0,0],
                        [0,0,1,0],[SxB2,0,0,F]])
    @property
    def A_b(self):
        '''F coordinate constraint on W momentum: ellipsoid'''
        K,B=self.K, self.b.Beta()
        mW2, x0p=self.mW2, self.x0p
        A_b_=np.array([1-B*B,0,0,B*x0p],[0,1,0,0],[0,0,1,0],[B*x0p,0,0,mW2-x0p**2])
        return K.dot(A_b_).dot(K.T)
    @property
    def R_T(self):
        '''Rotation from F coordinat to laboratory coordinate'''
        b_xyz=self.b.Px(),self.b.Py(),self.b.Pz()
        R_z=R(2,-self.mu.Phi())
        R_y=R(1,0.5*math.pi-self.mu.Theta())
        R_x=next(R(0,-math.atan2(z,y)) for x,y,z in (R_y.dot(R_z.dot(b_xyz)),))
        return R_z.T.dot(R_y.T.dot(R_x.T))
    @property
    def H_tilde(self):
        '''Transformation of t=[c,s,1] to p_nu: F coordinate'''
        x1, y1, p = self.x1,self.y1,self.mu.P()
        Z, w, Om = self.Z, self.w, math.sqrt(self.Om2)
        return np.array([[Z/Om,0,x1-p],[w*Z/Om,0,y1],[0,Z,0]])
    @property
    def H(self):
        '''Transformation of t=[c,s,1] to: lab coordinate.'''
        return self.R_T.dot(self.H_tilde)
    @property
    def H_perp(self):
        '''Transformation of t=[c,s,1] to pT_nu: lab coordinate.'''
        return np.vstack(self.H[:2],[0,0,1])
    @property
    def N(self):
        '''Solution ellipse of pT_nu: lab coordinate.'''
        HpInv=np.linalg.inv(self.H_perp)
        return HpInv.T.dot(UnitCircle()).dot(HpInv)

class singleNeutrinoSolution(object):
    '''Most likely neutrino momentum for tt-->lepton+jets'''
    def __init__(self, b , mu, #Lorentz vectors
               metX, metY,   #Momentum imbalance
               sigma2,          #Momentum imlance unc. 2X2 matrix
              mW2=mW**2,mT2=mT**2):
        self.solutionSet=nuSolutionSet(b,mu,mW2,mT2)
        S2=np.vstack([np.vstack([np.linalg.inv(sigma2),[0,0]]).T,[0,0,0]])
        V0 = np.outer([metX, metY, 0], [0, 0, 1])
        deltaNu=V0-self.solutionSet.H
        self.X=np.dot(deltaNu.T,S2).dot(deltaNu)
        M=next(XD+XD.T for XD in (self.X.dot(Derivative()),))
        solutions=intersections_ellipses(M,UnitCircle())
        self.solutions=sorted(solutions,key=self.calcX2)
    def calcX2(self,t):
        return np.dot(t,self.X).dot(t)
    @property
    def chi2(self):
        return self.calcX2(self.solutions[0])
    @property
    def nu(self):
        '''Solution for neutrino momentum'''
        return self.solutionSet.H.dot(self.solutions[0])
class doubleNeutrinoSolutions(object):
    '''Solution pairs of neutrino momenta, tt-->leptons'''
    def __init__(self, b,b_ , mu,mu_, #Lorentz vectors
               metX, metY,   #Momentum imbalance
               sigma2,          #Momentum imlance unc. 2X2 matrix
              mW2=mW**2,mT2=mT**2):
        self.solutionSet=[nuSolutionSet(B,M,mW2,mT2)
                          for B, M in zip((b,b_),(mu,mu_))]
        V0=np.outer([metX, metY, 0],[0,0,1])
        self.S=V0-UnitCircle()
        N,N_=[ss.N for ss in self.solutionSets]
        n_=self.S.T.dot(N_).dot(self.S)
        v=intersections_ellipses(N,n_)
        v_=[self.S.dot(sol) for sol in v]
        if not v and leastsq:
            es=[ss.H_perp for ss in self.solutionSets]
            met=np.array([metX,metY,1])
            def nus(ts):
                return tuple(e.dot([math.cos(t),math.sin(t),1])
                             for e, t in zip(es, ts))
            def residulas(params):
                return sum(nus(paras),-met)[:2]
            ts, _=leastsq(residuals,[0,0],ftol=5e-5,epsfcn=0.01)
            v,v_=[[i] for i in nus(ts)]
        for k, v in {'perp': v,"perp_":v_,"n_":n_}.items():
            setattr(self, k, v)
    @property
    def nunu_s(self):
        '''Solution pairs for neutrino momenta'''
        K, K_=[ss.H.dot(np.linalg.inv(ss.H_perp)) 
               for ss in self.solutionSets ]
        return [(K.dot(s),K_.dot(s_)) 
               for s, s_ in zip(self.perp,self.perp_)]
'''//////////////////////////////////////////////////////////////'''
'''recontruct top quark pairs in semileptonic channel'''
from ROOT import TFile,TTree, gROOT
from ROOT import TLorentzVector, Math
from math import sqrt,log, cos,sin
from array import array
#define a function to reconstruct top quark pairs
xi_thad=37
x0_thad=215
xi_wlep=3.4
x0_wlep=78
xi_tlep=8.4
x0_tlep=176
def recons_tt(jets_p4, jets_btag, #4-vector and btag of the leading four jets 
              lepton_p4, #4-vector of the lepton
              sigma2, #Momentum imlance unc. 2X2 matrix
              metX, metY):
    kk=0 #for bjets index
    tt=0 # for light jets index 
    bjets_index=2*[0]
    jets_index=2*[0]
    for ijet in range(4):
        if(jets_btag[ijet]==1 and kk <2):
              bjets_index[kk]=ijet
              kk=kk+1
        else:
              jets_index[tt]=ijet
              tt=tt+1
    likelihoods=[]
    neutrino_p4=[]
    for i in range(2):
        bl=bjets_index[i]
        bh=bjets_index[i+1 if i==0 else i-1]
        singleNuSol=singleNeutrinoSolution(jets_p4[bl],lepton_p4,metX,metY,
                                           sigma2,mW**2,mT**2)
        nu_p3=singleNuSol.nu
        nu_E=sqrt(nu_p3[0]**2+nu_p3[1]**2+nu_p3[2]**2)
        nu_p4=TLorentzVector(nu_p3[0],nu_p3[1],nu_p3[2],nu_E)
        neutrino_p4.append(nu_p4)
        mass_wlep=(nu_p4+lepton_p4).M()
        mass_whad=(jets_p4[jets_index[0]]+jets_p4[jets_index[1]]).M()
        mass_tlep=(nu_p4+jets_p4[bl]+lepton_p4).M()
        mass_thad=(jets_p4[jets_index[0]]+jets_p4[jets_index[1]]+jets_p4[bh]).M()
        pro_wlep=Math.landau_pdf(mass_wlep,xi_wlep,x0_wlep)
        pro_tlep=Math.landau_pdf(mass_tlep,xi_tlep,x0_tlep)
        pro_thad=Math.landau_pdf(mass_thad,xi_thad,x0_thad)
        likelihood=-log(pro_wlep)-log(pro_tlep)-log(pro_thad)
        likelihoods.append(likelihood)
    if (likelihoods[0] < likelihoods[1]):
        bl=bjets_index[0]
        nu_p4=neutrino_p4[0]
        bh=bjets_index[1]

    else:
        bl=bjets_index[1]
        nu_p4=neutrino_p4[1]
        bh=bjets_index[0]
    wl_p4=lepton_p4+nu_p4
    wh_p4=jets_p4[jets_index[0]]+jets_p4[jets_index[1]]
    topl_p4=jets_p4[bl]+lepton_p4+nu_p4
    toph_p4=jets_p4[bh]+jets_p4[jets_index[0]]+jets_p4[jets_index[1]]
    back_Dict={"topl_p4":topl_p4,"toph_p4":toph_p4,"nu_p4":nu_p4,
              "bl_index":bl,"bh_index":bh,"wh_p4":wh_p4,"wl_p4":wl_p4}
    return backDict
'''///////////////////////////////////////////////////////////////////'''
inputFile="new_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root"
#evetn selection has been applied in get_info.cpp
file=TFile(inputFile,"update")
mytree=file.Get("mytree") # get the original tree
Nevents = mytree.GetEntries()
jnum=4  #define the number of jets needed to reconstruct ttbar
sigma2=[[1.1,0.0],[0.0,1.1]]
for i in range(10):
    mytree.GetEntry(i)
    metX=mytree.MET_pt*cos(mytree.MET_phi);
    metY=mytree.MET_pt*sin(mytree.MET_phi);
    lep_eta=mytree.lepton_eta[0]
    lep_pt=mytree.lepton_pt[0]
    lep_mass=mytree.lepton_mass[0]
    lep_phi=mytree.lepton_phi[0]
    lep_p4=TLorentzVector() #4-vector of lepton
    lep_p4.SetPtEtaPhiM(lep_pt,lep_eta,lep_phi,lep_mass)
    mon_jets=[] #store 4-vector of the leading 4 jets
    jets_btag=[] # store btag(1 or 0) of the leading 4 jets
    
    for k in range(jnum):
        jet_p4=TLorentzVector()
        jet_pt=mytree.Jet_pt[k]
        jet_eta=mytree.Jet_eta[k]
        jet_phi=mytree.Jet_phi[k]
        jet_mass=mytree.Jet_mass[k]
        jet_p4.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass)
        mon_jets.append(jet_p4)
        jets_btag.append(mytree.Jet_btaged[k])
    #end the loop over the four leading jets
    #def recons_tt(jets_p4, jets_btag, lepton_p4, sigma2, metX, metY)
    
    nuslo=singleNeutrinoSolution(mon_jets[0],lep_p4,metX,metY,sigma2,mW**2,mT**2)
    nu_p3=nuslo.nu
    print(nu_p3[0],nu_p3[1],nu_p3[2])
    print("-----------------")

    '''
    tt_dict=recons_tt(mon_jets,jets_btag,lep_p4,sigma2,metX,metY)
    nu_p4=tt_dict["nu_p4"]      #4-vector of neutrino
    topl_p4=tt_dict["topl_p4"]  #4-vector of leptonic top
    toph_p4=tt_dict["toph_p4"]   #4-vector of hadronic top
    bl_index=tt_dict["bl_index"] #leptonic bjet index
    bh_index=tt_dict["bh_index"] #hadronic bjet index
    wl_p4=tt_dict["wl_p4"]
    wh_p4=tt_dict["wh_p4"]
    print(nu_p4.Px()," ", nu_p4.Py(),nu.Pz())
  '''











