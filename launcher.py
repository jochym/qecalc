from setting import Setting

class Launcher(Setting):
    def __init__(self, fname=None):
        Setting.__init__(self,fname)

    def check(self, x):
        """Will check the exit status of QE program"""
        signal = x & 0xFF
        exitcode = (x >> 8) & 0xFF
        if exitcode != 0:
            raise Exception("Quantum Espresso crashed: check your settings")

    def pwscfLauncher(self):
        import os
        cmdstr = self.paraPrefix + " pw.x " +  self.paraPostfix + " -inp " + \
                 self.pwscfInput + " > " + self.pwscfOutput + "< /dev/null"
        print cmdstr         
        self.check(os.system(cmdstr))

    def singlePhononLauncher(self):
        import os
        self.pwscfLauncher()
        cmdstr_ph = self.paraPrefix + " ph.x " +  self.paraPostfix + " -inp " + \
                 self.phInput + " > " + self.phOutput + "< /dev/null"
        print cmdstr_ph        
        self.check(os.system(cmdstr_ph))

        cmdstr_dynmat = "dynmat.x < " + self.dynmatInput + " > " + self.dynmatOutput
        print cmdstr_dynmat
        self.check(os.system(cmdstr_dynmat))
        
    def multiPhononLauncher(self):
        import os
        self.pwscfLauncher()
        cmdstr_ph = self.paraPrefix + " ph.x " +  self.paraPostfix + " -inp " + \
                 self.phInput + " > " + self.phOutput + "< /dev/null"
        print cmdstr_ph        
        self.check(os.system(cmdstr_ph))
        cmdstr_q2r = "q2r.x < " + self.q2rInput + " > " + self.q2rOutput
        print cmdstr_q2r
        self.check(os.system(cmdstr_q2r))
        cmdstr_matdyn = "matdyn.x < " + self.matdynInput + " > " + self.matdynOutput
        print cmdstr_matdyn
        self.check(os.system(cmdstr_matdyn))
