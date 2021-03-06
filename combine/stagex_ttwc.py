from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import fnmatch 

class stagex_ttwc(PhysicsModel):
    "Allow different signal strength fits for the stage-x model"
    def __init__(self):
		self.POIs = ""
		self.acttbar=False 
		self.acggh=False 
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
        	print "this is ttbar process"
		muname=""
		if self.stage0:
			if "ttbar" in process:
				if self.acttbar:
					if "ci0001" in process:
						muname="r_tt_times_ci0001"
					elif "ci2000" in process:
						muname="r_tt_times_ci2000"
					elif "ci1100" in process:
						muname="r_tt_times_ci1100"
					elif "ci1000" in process:
						muname="r_tt_times_ci1000"
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
			self.modelBuilder.doVar("x[0.0,-5,5]" )
			self.modelBuilder.doVar("y[0.0,-5,5]" )
			self.modelBuilder.doVar("z[0.0,-5,5]" )
			self.modelBuilder.doVar("k[0.0,-5,5]" )
			# x:Cpq3, y: Cpu, z:ReCup, k: ImCup
			self.modelBuilder.factory_("expr::r_tt_times_ci0001(\"@0*@0\", k)")
			self.modelBuilder.factory_("expr::r_tt_times_ci2000(\"-0.5*@0+0.5*@0*@0\", x)")
			self.modelBuilder.factory_("expr::r_tt_times_ci1100(\"@0*@1\", x,y)")
			self.modelBuilder.factory_("expr::r_tt_times_ci1000(\"2*@0-@0*@1-@0*@0\", x, y)")
			self.modelBuilder.factory_("expr::r_tt_times_ci0100(\"2*@0-@0*@1-@0*@0\", y,x)")
			self.modelBuilder.factory_("expr::r_tt_times_ci0200(\"-0.5*@0+0.5*@0*@0\", y)")
			self.modelBuilder.factory_("expr::r_tt_times_ci0010(\"2*@0-@0*@0\", z)")
			self.modelBuilder.factory_("expr::r_tt_times_ci0020(\"-0.5*@0+0.5*@0*@0\", z)")
			self.modelBuilder.factory_("expr::r_tt_times_ci0000(\"1-@0*@0-1.5*@1+0.5*@1*@1-1.5*@2+@1*@2+0.5*@2*@2-1.5*@3+0.5*@3*@3\", k,x,y,z)")

			self.pois.append("x,y,z,k")
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
		    if 'doacggh' in po: 
			    self.acggh= True
			    print "doing tth AC ggH"
		    if 'singlemu' in po: 
			    self.singlemu= True
			    print "doing single mu"

    def doParametersOfInterest(self):
	    self.POIs=",".join(self.pois)
	    print "the parameters of interest: ", self.POIs
	    self.modelBuilder.doSet("POI",self.POIs)

stagex_ttwc = stagex_ttwc()
