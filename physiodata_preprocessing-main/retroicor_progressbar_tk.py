from numpy import *
import pylab
import statsmodels.api as sm
#import nifti
import nibabel as nib
import os
import time

from Tkinter import *
import Tkconstants, tkFileDialog

from js_validateentry import *
#from progressmeter import *

import progressbar

class App:
    def __init__(self,master,title = 'Retroicor') :
        self.frame = Frame(master)
        self.frame.pack()
        self.TR=2.0
        self.SampFreq=200
        self.AdditionalDelay = 0.0

        images=Button(self.frame, text="EPIFile",  width=18,  command = self.EPIFile).grid()

        resp = Button(self.frame, text= 'Respiration File', width=18,  command = self.RespFile).grid(row=1)

        card = Button(self.frame, text= 'Pulse File', width=18,  command = self.CardioFile).grid(row=2)

        go = Button(self.frame,  text='GO', command = self.start).grid(row=100)
        cn = Button(self.frame, text = 'Exit', command = root.destroy).grid(row=100, column = 1)
        #self.status=Meter(self.frame)
        #self.status.grid(row=200, columnspan=2)

        label_tr=Label(self.frame, text='TR [s]', width = 22,  anchor='e')
        label_tr.grid(row=3)
        self.TR = FloatEntry(self.frame, str(self.TR))
        self.TR.grid(row=3, column=1)

        label_samp=Label(self.frame, text='Sampling Freq. Physio [Hz]', anchor='e',  width = 22)
        label_samp.grid(row=4)
        self.SampFreq = IntegerEntry(self.frame, str(self.SampFreq))
        self.SampFreq.grid(row=4, column=1)

        label_delay=Label(self.frame, text='Additional Delay cp only [s]', anchor='e',  width = 22)
        label_delay.grid(row=5)
        self.AdditionalDelay = FloatEntry(self.frame, str(self.AdditionalDelay))
        self.AdditionalDelay.grid(row=5, column=1)

        Label(self.frame, text='Acquisition order', width=22,  anchor='e').grid(row=6)
        self.listbox = Listbox(self.frame, height = 4)
        for item in ["Ascending", "Descending", "Interleaved", "Interleaved 2"]:
            self.listbox.insert(END, item)
        self.listbox.grid(row=6, column=1)


    def EPIFile(self):
        '''Load MR data Nifti Format'''
        self.MRfilename = tkFileDialog.askopenfilename(filetypes=[('all files','*.*'),('NIFTI','*.nii'), ('NIFTI.GZ', '*.nii.gz')], title='EPI File')
        if self.MRfilename:
            #try: self.MRdata=nifti.NiftiImage(self.MRfilename)
            try: self.MRdata=nib.load(self.MRfilename)
            except IOError: sys.exit("File not found") 
            self.hdr=self.MRdata.header
            #self.fdata=self.MRdata.get_fdata()
            #self.Dims=self.fdata.shape # t,z,y,x
            self.Dims=self.MRdata.shape # x,y,z,t
            Label(self.frame,  text = os.path.basename(self.MRfilename) + ': ' + str(self.Dims[3]) + 'Vols' ).grid(row=0, column=1)


    def RespFile(self):
        '''Load respiratory phase data from tk_physionoise.py'''
        self.Respfilename = tkFileDialog.askopenfilename(filetypes=[('txt files','*.txt')], title='Respiration File')
        if self.Respfilename:
            try: self.Resp=fromfile(self.Respfilename,'f', sep=' ' )
            except IOError: sys.exit("Respiration File not found")
            Label(self.frame,  text = os.path.basename(self.Respfilename)).grid(row=1, column=1)


    def CardioFile(self):
        '''Load Pulse phase data from tk_physionoise.py'''
        self.Cardiofilename = tkFileDialog.askopenfilename(filetypes=[('txt files','*.txt')], title='Cardio File')
        if self.Cardiofilename:
            try: self.Cardio=fromfile(self.Cardiofilename,'f', sep=' ' )
            except IOError: sys.exit("Cardiac File not found")
            Label(self.frame,  text = os.path.basename(self.Cardiofilename)).grid(row=2, column=1)


    def check(self):
        '''Check ob alle Files vorhanden sind und Acq Order ausgewaehlt ist'''
        try:
            if self.MRfilename and self.Respfilename and self.Cardiofilename and self.listbox.curselection():
                return True
        except AttributeError:
            return False

    def sliceorder(self):
        selection = self.listbox.get(self.listbox.curselection())
        if selection == 'Ascending':
            self.Sorder = range(self.Dims[2])
        elif selection == 'Descending':
            self.Sorder = range(self.Dims[2]-1, -1, -1)
        elif selection == 'Interleaved':
            self.Sorder = range(1, self.Dims[2], 2) + range(0, self.Dims[2], 2)
        elif selection == 'Interleaved 2':
            self.Sorder = range(0, self.Dims[2], 2) + range(1, self.Dims[2], 2)


    def start(self):
        if self.check():
            TR=float(self.TR.results.get())
            SampFreq = int(self.SampFreq.results.get())
            AdditionalDelay = float(self.AdditionalDelay.results.get())
            SliceTR= TR / float(self.Dims[2])  # [s] Acqu. time per slice
            self.sliceorder()
            slices=self.MRdata.dataobj

            rp=zeros((self.Dims[3], self.Dims[2]))
            cp=zeros((self.Dims[3], self.Dims[2]))

            for sl in range(self.Dims[2]):
                rp[:,self.Sorder[sl]]=self.Resp[int(rint((((sl + 0.5)*SliceTR) ) * SampFreq)): -1: int(TR*SampFreq)][0:self.Dims[3]]
                cp[:,self.Sorder[sl]]=self.Cardio[int(rint((((sl + 0.5)*SliceTR) + AdditionalDelay) * SampFreq)): -1: int(TR*SampFreq)][0:self.Dims[3]]

            residual=zeros((self.Dims))
            start = time.clock()
            with progressbar.ProgressBar(self.Dims[2]) as bar:
                for sl in range(self.Dims[2]):
                    time.sleep(0.1) #progressbar
                    bar.update(value = (sl/float(self.Dims[2]-1))) #progressbar
                    X=sm.tools.add_constant(column_stack((sin(rp[:, sl]), cos(rp[:, sl]) , sin(2*rp[:, sl]),  cos(2*rp[:, sl]),  sin(cp[:, sl]),  cos(cp[:, sl]), sin(2*cp[:, sl]),  cos(2*cp[:, sl]))))
                    for y in range(self.Dims[1]):
                        for x in range(self.Dims[0]):
                            erg = sm.GLS(slices[x, y, sl, :], X).fit()
                            residual[x, y, sl, :] = erg.resid + erg.params[0] #self.status.set(value = (sl/float(self.Dims[2]-1)) ) # progressmeter
                            self.frame.update()
                            #print(erg.summary())
                elapsed = (time.clock() - start)
                print "Dauer :" + str(elapsed)    
            #Save residuals
            outdata=nib.Nifti1Image(residual, affine=self.MRdata.affine, header=self.hdr)
            #outdata.header=self.MRdata.header
            # Dateiendung abschneiden .nii oder nii.gz
            name, ext=os.path.splitext(self.MRfilename)
            # wenn gz dann auch nii abschneiden
            if (ext=='.gz'):
                name, ext=os.path.splitext(name)
            outfilename = name + '_retroicor_' + str(AdditionalDelay * 1000) + '_ms.nii.gz'
            #outdata.save(outfilename)
            outdata.to_filename(os.path.join('build',outfilename))



if __name__=='__main__':

    root = Tk()
    root.title('Retroicor')
    app=App(root)
    root.mainloop()
