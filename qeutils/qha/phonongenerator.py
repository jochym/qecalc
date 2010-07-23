import os
from qecalc.multiphononcalc import MultiPhononCalc

class PhononGenerator:
    def __init__(self, multiPhononCalc, pickleName = None):
        
       self._mphon = multiPhononCalc
       if pickleName == None:
           pickleName = 'volume_expansion.pkl'
       self.pickleName = pickleName
    
       
    def load(self, pickleName = None, input = None):
        """
        Loads output from VolumeOptimizer either as a pickle or a dictionary 
        """
        if input != None:
            self._input = input
        else:            
            import pickle
            if pickleName == None:
                pickleName = self.pickleName
            self._input = pickle.load( open(pickleName, 'r') )
        
        
    def volume(self):
        """
        Volume expansion in percents
        """
        return self._input["volume"]
      
       
    def launch(self, pickleName = None, input = None, volume = None, prefix = ''):
        """
        Performs series of multiple phonon calculations at provided volume expansions
        Loads output from VolumeOptimizer either as a pickle or a dictionary (input)
        volume - list of volume expansions (default - all available)        
        """
        
        self.load(pickleName, input)
        if volume == None:
            volume = self.volume()
        
        for v in  volume:
            if v in self.volume():
                i = self.volume().index(v)
                self._mphon.pw.input.parse()
                self._mphon.pw.input.structure.load( source = 'pwconfig', \
                                               configString = self._config()[i] )
                
                #print self._mphon.pw.input.toString()
                self._mphon.pw.input.save()
                self._mphon.launch()
                self._pack( prefix = str(i) + prefix)
                
                                
    def _config(self):
        return self._input["config"]
    
    
    def _pack(self, prefix):
        fileList = []
        for task in self._mphon.taskList:
            task.syncSetting()
            fileList.extend( task.setting.getExistingFiles() )
        
        #solve for duplicates:
        fileDic = {}
        for name in fileList:
            fileDic[name] = name
        fileStr = ''
        for name in fileDic:
            fileStr = fileStr + ' ' + str(name)
        
        print fileStr
        
        archiveName = str(prefix) + '.tgz'
        os.system('tar -zcf ' + archiveName + ' ' + fileStr )
        
        
if __name__ == '__main__':

    PhononGenerator(multiPhononCalc = MultiPhononCalc(filename = 'config.ini')).launch()
        
