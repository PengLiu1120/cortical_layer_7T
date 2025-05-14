from numpy import *
import pylab
#import scikits.statsmodels as sm
import statsmodels.api as sm
#import nifti
import nibabel as nib
import os
import time

from Tkinter import *
import Tkconstants, tkFileDialog

from js_validateentry import *
#from progressmeter import *

class App:
	def __init__(self,master,title = 'Retroicor') :
		self.frame = Frame(master)
		self.frame.pack()
		self.TR=2.0
		self.SampFreq=200
		self.AdditionalDelay = 0.025
		self.timepoints = 1
		self.step = 0.0

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

		label_delay=Label(self.frame, text='Additional Delay [s]', anchor='e',  width = 22)
		label_delay.grid(row=5)
		self.AdditionalDelay = FloatEntry(self.frame, str(self.AdditionalDelay))
		self.AdditionalDelay.grid(row=5, column=1)

		Label(self.frame, text='Acquisition order', width=22,  anchor='e').grid(row=6)
		self.listbox = Listbox(self.frame, height = 4)
		for item in ["Ascending", "Descending", "Interleaved", "Interleaved 2"]:
			self.listbox.insert(END, item)
		self.listbox.grid(row=6, column=1)

		Label(self.frame,  text="Timepoints",  width = 22,  anchor='e').grid(row=10)
		self.timepoints=IntegerEntry(self.frame, str(self.timepoints))
		self.timepoints.grid(row=10,  column=1)
        
		Label(self.frame,  text='Step [s]',  width = 22,  anchor='e').grid(row=11)
		self.step=FloatEntry(self.frame,  str(self.step))
		self.step.grid(row=11,  column=1)

	def EPIFile(self):
		'''Load MR data Nifti Format'''
		self.MRfilename = tkFileDialog.askopenfilename(filetypes=[('all files','*.*'),('NIFTI','*.nii'), ('NIFTI.GZ', '*.nii.gz')], title='EPI File')
		if self.MRfilename:
			#try: self.MRdata=nifti.NiftiImage(self.MRfilename)
			try: self.MRdata=nib.load(self.MRfilename)
			except IOError: sys.exit("File not found")
			self.hdr=nib.Nifti1Header() 
			#self.Dims=self.MRdata.data.shape # t,z,y,x
			self.Dims=self.MRdata.shape # x,y,z,t
			#Label(self.frame,  text = os.path.basename(self.MRfilename) + ': ' + str(self.Dims[0]) + 'Vols' ).grid(row=0, column=1)
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
	# def sliceorder(self):
	#     selection = self.listbox.get(self.listbox.curselection())
	#     if selection == 'Ascending':
	#         self.Sorder = range(self.Dims[1])
	#     elif selection == 'Descending':
	#         self.Sorder = range(self.Dims[1]-1, -1, -1)
	#     elif selection == 'Interleaved':
	#         self.Sorder = range(1, self.Dims[1], 2) + range(0, self.Dims[1], 2)
	#     elif selection == 'Interleaved 2':
	#         self.Sorder = range(0, self.Dims[1], 2) + range(1, self.Dims[1], 2)
		
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
			#SliceTR= TR / self.Dims[1]  # [s] Acqu. time per slice
			SliceTR= TR / float(self.Dims[2])  # [s] Acqu. time per slice
			timepoints = int(self.timepoints.results.get())
			step = float(self.step.results.get())
            
			self.sliceorder()
			#slices=self.MRdata.data
			slices=self.MRdata.dataobj

			#rp=zeros((timepoints, self.Dims[0], self.Dims[1]))
			#cp=zeros((timepoints, self.Dims[0], self.Dims[1]))

			rp=zeros((timepoints, self.Dims[3], self.Dims[2]))
			cp=zeros((timepoints, self.Dims[3], self.Dims[2]))
            
			for delay in range(timepoints):
				for sl in range(self.Dims[2]):
					rp[delay, :, self.Sorder[sl]]=self.Resp[int((((sl + 0.5)*SliceTR) ) * SampFreq): -1: int(TR*SampFreq)][0:self.Dims[3]]
					cp[delay, :, self.Sorder[sl]]=self.Cardio[int((((sl + 0.5)*SliceTR) + delay * step + AdditionalDelay) * SampFreq): -1: int(TR*SampFreq)][0:self.Dims[3]]

			#residual=zeros((timepoints, self.Dims[1], self.Dims[2], self.Dims[3]))
			#residual=zeros((timepoints, self.Dims[2], self.Dims[1], self.Dims[0]))
			residual=zeros((self.Dims[0], self.Dims[1], self.Dims[2], timepoints))

			for delay in range(timepoints):
				#for sl in range(self.Dims[1]):
				for sl in range(self.Dims[2]):	
					#X=sm.tools.add_constant(cp[delay, :, sl])
					#X=sm.tools.add_constant(column_stack((sin(rp[delay, :, sl]), cos(rp[delay, :, sl]) , sin(2*rp[delay, :, sl]),  cos(2*rp[delay, :, sl]),  sin(cp[delay, :, sl]),  cos(cp[delay, :, sl]), sin(2*cp[delay, :, sl]),  cos(2*cp[delay, :, sl]))))
					X=sm.tools.add_constant(column_stack((sin(cp[delay, :, sl]),  cos(cp[delay, :, sl]), sin(2*cp[delay, :, sl]),  cos(2*cp[delay, :, sl]))))
					#for y in range(self.Dims[2]):
						#for x in range(self.Dims[3]):
					for y in range(self.Dims[1]):
						for x in range(self.Dims[0]):		                     	
							#erg = sm.GLS(slices[:, sl, y, x], X).fit()
							erg = sm.GLS(slices[x, y, sl, :], X).fit()
							#residual[delay, sl, y, x] = erg.fvalue	#correct?
							residual[x, y, sl, delay] = erg.fvalue
				#self.status.set(value = ((delay +1.) / timepoints  ) ) # progressmeter
				self.frame.update()
			#Save residuals
			#outdata=nifti.NiftiImage(residual)
			self.hdr['dim'] = [4, 104, 104, timepoints, 1, 1, 1, 1]
			outdata=nib.Nifti1Image(residual, affine=self.MRdata.affine, header=self.hdr)
			#outdata.header=self.hdr
			#outdata.header['dim'] = [4, 106, 106, timepoints, 1, 1, 1, 1]
			# Dateiendung abschneiden .nii oder nii.gz
			name, ext=os.path.splitext(self.MRfilename)
			# wenn gz dann auch nii abschneiden
			if (ext=='.gz'):
				name, ext=os.path.splitext(name)
			#outfilename = name + '_retroicor_' + str(AdditionalDelay * 1000) + '_timepoints_cp_only.nii.gz'
			outfilename = name + '_retroicor_' +  '20steps_50ms_cardio_only.nii.gz'
			#outdata.save(outfilename)
			outdata.to_filename(os.path.join('build',outfilename))

			#nib.save(outdata,outfilename)


if __name__=='__main__':

    root = Tk()
    root.title('Retroicor')
    app=App(root)
    root.mainloop()
