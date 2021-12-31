from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import fnmatch
from math import sqrt,copysign 
sign = lambda x: copysign(1, x) 

class stagex_ttwc3(PhysicsModel):
    "Allow different signal strength fits for the stage-x model"
    def __init__(self):
		self.POIs = ""
		self.acttbar=False 
		self.acfcp=False 
		self.accor=False 
		self.options= ""
		self.stage0= False 
		self.rvrf= False 
		self.singlemu= False 
		self.splitHad= False 
		self.muNames =[]
		self.pois = []
		self.count=0
    def setModelBuilder(self, modelBuilder):
		PhysicsModel.setModelBuilder(self, modelBuilder)
		self.modelBuilder.doModelBOnly = False
    def getYieldScale(self,bin,process):
        if "ttbar" not in process: # return 1 for background
		    return 1
        else:
        	print process, "is a ttbar process with EW correction"
		if self.stage0:
			if "ttbar" in process:
				if self.acttbar:
					if "ci0001" in process:
						muname="r_tt_times_ci0001"
					elif "ci0100" in process:
						muname="r_tt_times_ci0100"
					elif "ci0200" in process:
						muname="r_tt_times_ci0200"
					elif "ci0010" in process:
						muname="r_tt_times_ci0010"
					elif "ci0020" in process:
						muname="r_tt_times_ci0020"
					elif "ci0000" in process:
						muname="r_tt_times_ci0000"
				else:
					muname="r_tt"
				print "muname: %s" %muname
	
        if self.acttbar:
			if not self.acfcp:
				# x:Cpq3, y: Cpu, z:ReCup, k: ImCup
				self.modelBuilder.doVar("y[0.0,-3,3]" )
				self.modelBuilder.doVar("z[0.0,-3,3]" )
				self.modelBuilder.doVar("k[0.0,-3,3]" )
				self.modelBuilder.factory_("expr::r_tt_times_ci0001(\"@0*@0\", k)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0100(\"2*@0-@0*@0\", y)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0200(\"-0.5*@0+0.5*@0*@0\", y)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0010(\"2*@0-@0*@0\", z)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0020(\"-0.5*@0+0.5*@0*@0\", z)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0000(\"1-@0*@0-1.5*@1+0.5*@1*@1-1.5*@2+0.5*@2*@2\", k,y,z)")
				self.pois.append("y,z,k")
			else:
				self.modelBuilder.doVar("y[0.0,-3,3]" )
				self.modelBuilder.doVar("fcp[0.0,-1.0,1.0]")
				self.modelBuilder.doVar("mu[1.0,-10,10]")
				self.modelBuilder.factory_("expr::r_tt_times_ci0001(\"pow(sign(@0)*sign(@1)*sqrt(abs(@0)*abs(@1)),2)\", fcp,mu)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0100(\"2*@0-@0*@0\", y)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0200(\"-0.5*@0+0.5*@0*@0\", y)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0010(\"2*(1-sign(@1)*sqrt((1-abs(@0))*abs(@1)))\
-pow(1-sign(@1)*sqrt((1-abs(@0))*abs(@1)),2)\", fcp,mu)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0020(\"-0.5*(1-sign(@1)*sqrt((1-abs(@0))*abs(@1)))+\
0.5*pow((1-sign(@1)*sqrt((1-abs(@0))*abs(@1))),2)\", fcp,mu)")
				self.modelBuilder.factory_("expr::r_tt_times_ci0000(\"1-pow(sign(@0)*sign(@1)*sqrt(abs(@0)*abs(@1)),2)\
-1.5*@2+0.5*@2*@2-1.5*(1-sign(@1)*sqrt((1-abs(@0))*abs(@1)))\
+0.5*pow((1-sign(@1)*sqrt((1-abs(@0))*abs(@1))),2)\", fcp,mu,y)")
				self.pois.append("y,fcp,mu")
			self.POIs=",".join(self.pois)
			self.modelBuilder.doSet("POI",self.POIs)
			print "parameters of interest in ttbar analysis: ", self.POIs
        	
        if self.modelBuilder.out.var(muname):
			print "reclying %s" %muname
        else:                              
			if not "times" in muname:
				self.modelBuilder.doVar("%s[1,0,10]" % muname)
				print "scale process %s with %s" %(process,muname)
				self.pois.append(muname)
				self.POIs=",".join(self.pois)
				self.modelBuilder.doSet("POI",self.POIs)
				print "Default parameters of interest: ", self.POIs
				
        return muname 
    def setPhysicsOptions(self,physOptions):
	    for po in physOptions:
		    if 'doStage0' in po: 
			    self.stage0= True
			    print "doing stage0"
		    if 'doacttbar' in po: 
			    self.acttbar= True
			    print "doing tth AC ttbar"
		    if 'dofcp' in po: 
			    self.acfcp= True
			    print "doing tth AC ggH"
		    if 'singlemu' in po: 
			    self.singlemu= True
			    print "doing single mu"

    def doParametersOfInterest(self):
	    self.POIs=",".join(self.pois)
	    print "the parameters of interest: ", self.POIs
	    self.modelBuilder.doSet("POI",self.POIs)

stagex_ttwc3 = stagex_ttwc3()
