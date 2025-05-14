#!/home/juli/python-test/bin/python
# -*- coding: utf-8 -*-
import os 
import subprocess
from Tkinter import *
import Tkconstants, tkFileDialog

programm='/home/doehlerj/Dokumente/PhD/sensemap/fMRI/physio/JS_PhysioNoise_pyLog.py'

print(os.listdir(os.getcwd()))

#
# http://tkinter.unpythonic.net/wiki/ValidateEntry
#
class ValidatingEntry(Entry):
# base class for validating entry widgets

    def __init__(self, master, value="", **kw):
        apply(Entry.__init__, (self, master), kw)
        self.__value = value
        self.__variable = StringVar()
        self.__variable.set(value)
        self.__variable.trace("w", self.__callback)
        self.config(textvariable=self.__variable)
        self.results = StringVar()
        if self.__value is None: self.results.set(None)
        else:
            self.results.set(self.__value)
    
    def __callback(self, *dummy):
        value = self.__variable.get()
        newvalue = self.validate(value)
        if newvalue is None:
            self.__variable.set(self.__value)
        elif newvalue != value:
            self.__value = newvalue
            self.__variable.set(newvalue)
        else:
            self.__value = value
    
    def validate(self, value):
        # override: return value, new value, or None if invalid
        self.results.set(value)
        return value
    
    def getresults(self, value):
        # override: return value, or chopped value in the case of ChopLengthEntry
        return self.results.get()


'''
The first two examples are subclasses that check that the input is a valid Python integer or float, respectively.
The validate method simply tries to convert the value to an object of the right kind, and returns None (reject) if that fails.
'''

class IntegerEntry(ValidatingEntry):
    def validate(self, value):
        try:
            if value:
                v = int(value)
                self.results.set(value)
            return value
        except ValueError:
            return None


class FloatEntry(ValidatingEntry):
    def validate(self, value):
        try:
            if value:
                v = float(value)
                self.results.set(value)
            return value
        except ValueError:
            return None

class StringEntry(ValidatingEntry):
    #same as ValidatingEntry; nothing extra
    pass


class App:
    def __init__(self,master,title = 'PhysioNoise') :
        self.logfilename= tkFileDialog.askopenfilenames(filetypes=[('all files','*.*'),('log','*.txt')], title='Log File')
        
        Label(master, text='TR [s]').grid(row=0)
        self.TR = FloatEntry(master, str(2.0))
        self.TR.grid(row=0, column=1)

        #Label(master, text='output path').grid(row=0)
        #self.outputpath = StringEntry(master, str("/home/juli/Dokumente/PhD/sensemap/fMRI/physio/data/"))
        #self.outputpath.grid(row=0, column=1)
        
        Label(master, text='number of Volums').grid(row=0,  column=2)
        self.numTR=IntegerEntry(master,  str(130))
        self.numTR.grid(row=0, column=3)
        
        Label(master, text='Input Hz').grid(row=1,  column=0)
        self.inputHz = IntegerEntry(master,  str(200))
        self.inputHz.grid(row=1, column=1)
        
        Label(master, text='Highpass').grid(row=1,  column=2)
        self.fpass= IntegerEntry(master,  str(1))
        self.fpass.grid(row=1, column=3)
        
        Label(master, text='Lowpass').grid(row=2,  column=0)
        self.fstop=IntegerEntry(master,  str(10))
        self.fstop.grid(row=2,  column=1)
        
        Label(master, text='Trim Window').grid(row=2,  column=2)
        self.trimwindow = IntegerEntry(master,  str(20))
        self.trimwindow.grid(row=2, column=3)
        
        Label(master, text='plr').grid(row=3,  column=0)
        self.plr = FloatEntry(master,  str(1.2))
        self.plr.grid(row=3,  column=1)
        
        Label(master, text='StatThresh').grid(row=3,  column=2)
        self.statthresh = IntegerEntry(master,  str(6))
        self.statthresh.grid(row=3,  column=3)
        
        Label(master, text='rtimethresh').grid(row=4,  column=0)
        self.rtimethresh = IntegerEntry(master,  str(0))
        self.rtimethresh.grid(row=4, column=1)
        
        Label(master, text='cmagthresh').grid(row=4,  column=2)
        self.cmagthresh = IntegerEntry(master,  str(20))
        self.cmagthresh.grid(row=4,  column=3)
        
        Label(master, text='ctimethresh').grid(row=5,  column=0)
        self.ctimethresh = IntegerEntry(master,  str(0))
        self.ctimethresh.grid(row=5,  column=1)
        
#        Label(master, text='Prefix').grid(row=6,  column=0)
#        self.prefix=StringEntry(master,  'PreIcor')
#        self.prefix.grid(row=6, column=1)
        
        Label(master, text='Plots').grid(row=7)
        self.guiPlot=IntVar()
        Checkbutton(master, variable=self.guiPlot).grid(row=7, column=1)
          
        Label(master, text='Speed').grid(row=7,  column=2)
        self.guiSpeed=IntVar()
        Checkbutton(master, variable=self.guiSpeed).grid(row=7, column=3)
         
        ok= Button(master,  text ='Start',  command = self.run).grid(row=100)
        cn = Button(master, text = 'cancel', command = master.destroy).grid(row=100, column=1)
        
        
        

    def run(self):
        for fname in self.logfilename:
                parameter = ' -l ' + fname + ' --TR ' +  self.TR.results.get() + ' --numTR ' + self.numTR.results.get()\
                + ' --inputHz ' + self.inputHz.results.get() + ' --fpass ' + self.fpass.results.get() + ' --fstop ' + self.fstop.results.get()\
                + ' --trimwindow ' + self.trimwindow.results.get() + ' --plr ' + self.plr.results.get() + ' --statthresh ' \
                + self.statthresh.results.get() + ' --rtimethresh ' + self.rtimethresh.results.get() + ' --cmagthresh ' \
                + self.cmagthresh.results.get() + ' --ctimethresh ' + self.ctimethresh.results.get() + ' --prefix ' \
                + os.path.splitext(fname)[0]
        
                if self.guiPlot.get():
                    parameter +=' --plots '
                if self.guiSpeed.get():
                    parameter += ' --speed '
                
                print parameter
                os.chdir('/home/doehlerj/Dokumente/PhD/sensemap/fMRI/physio/data')
                print(os.getcwd())
                # os.system ( programm + parameter)
                subprocess.call('python ' + programm + parameter, shell=True)    
                print "Done"
if __name__=='__main__':
    root = Tk()
    app=App(root)
    root.mainloop()
