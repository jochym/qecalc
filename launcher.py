from setting import Setting

class Launcher(Setting):
    def __init__(self, fname=None):
        Setting.__init__(self,fname)

    def __check(self, x):
        """Will check the exit status of QE program"""
        signal = x & 0xFF
        exitcode = (x >> 8) & 0xFF
        if exitcode != 0:
            raise Exception("Quantum Espresso crashed: check your settings and/or clean your 'outdir' directory")

    def pwscfLauncher(self):
        import os
        cmdstr = self.paraPrefix + " pw.x " +  self.paraPostfix + " -inp " + \
                 self.pwscfInput + " > " + self.pwscfOutput + "< /dev/null"
        print cmdstr         
        self.__check(os.system(cmdstr))

    def phLauncher(self):
        import os
        cmdstr_ph = self.paraPrefix + " ph.x " +  self.paraPostfix + " -inp " + \
                 self.phInput + " > " + self.phOutput + "< /dev/null"
        print cmdstr_ph        
        self.__check(os.system(cmdstr_ph))

    def singlePhononLauncher(self):
        import os
        self.pwscfLauncher()
        self.phLauncher()

        cmdstr_dynmat = "dynmat.x -inp " + self.dynmatInput + " > " + self.dynmatOutput
        print cmdstr_dynmat
        self.__check(os.system(cmdstr_dynmat))

    def matdynLauncher(self):
        """Execute matdyn.x after successful run of pw.x + ph.x + q2r.x"""
        import os
        cmdstr_matdyn = "matdyn.x -inp " + self.matdynInput + " > " + self.matdynOutput
        print cmdstr_matdyn
        self.__check(os.system(cmdstr_matdyn))

    def multiPhononLauncher(self):
        """Runs complete sequence of programms needed to extract phonons except
        the last step: matdyn.x. Usecase: One then can regenerate matdyn.in
        for dispersions along different directions, phonon DOS etc"""
        import os
        self.pwscfLauncher()
        self.phLauncher()
        cmdstr_q2r = "q2r.x < " + self.q2rInput + " > " + self.q2rOutput
        print cmdstr_q2r
        self.__check(os.system(cmdstr_q2r))


    def multiPhononTaskLauncher(self):
        """Runs complete sequence of programms needed to extract phonons"""
        import os
        self.pwscfLauncher()
        self.phLauncher()
        cmdstr_q2r = "q2r.x < " + self.q2rInput + " > " + self.q2rOutput
        print cmdstr_q2r
        self.__check(os.system(cmdstr_q2r))
        cmdstr_matdyn = "matdyn.x < " + self.matdynInput + " > " + self.matdynOutput
        print cmdstr_matdyn
        self.__check(os.system(cmdstr_matdyn))
