import h5py
import numpy as np
from astropy.table import Table
from Payne.fitting import genmod

class getmod(object):
	"""docstring for genmod"""
	def __init__(self):
		super(getmod, self).__init__()

		self.photbands = ([
			'2MASS_H',
			'2MASS_J',
			'2MASS_Ks',
			'Bessell_B',
			'Bessell_I',
			'Bessell_R',
			'Bessell_U',
			'Bessell_V',
			'DECam_Y',
			'DECam_g',
			'DECam_i',
			'DECam_r',
			'DECam_u',
			'DECam_z',
			'GALEX_FUV',
			'GALEX_NUV',
			'Gaia_BP_DR2Rev',
			'Gaia_G_DR2Rev',
			'Gaia_RP_DR2Rev',
			'Hipparcos_Hp',
			'Kepler_D51',
			'Kepler_Kp',
			'PS_g',
			'PS_i',
			'PS_open',
			'PS_r',
			'PS_w',
			'PS_y',
			'PS_z',
			'SDSS_g',
			'SDSS_i',
			'SDSS_r',
			'SDSS_u',
			'SDSS_z',
			'TESS',
			'Tycho_B',
			'Tycho_V',
			'UKIDSS_H',
			'UKIDSS_J',
			'UKIDSS_K',
			'UKIDSS_Y',
			'UKIDSS_Z',
			'WISE_W1',
			'WISE_W2',
			'WISE_W3',
			'WISE_W4',
			])

		self.GM = genmod.GenMod()
		self.GM._initspecnn(nnpath='/Users/pcargile/Astro/ThePayne/Hecto_C3K_v2/C3K_Hecto_v2.h5')
		self.GM._initphotnn(self.photbands,nnpath='/Users/pcargile/Astro/MIST/nnMIST/')

	def gmod(self,*args,**kwargs):
		
		Teff = kwargs.get('Teff',5770.0)
		logg = kwargs.get('logg',4.44)
		FeH  = kwargs.get('FeH',0.0)
		aFe  = kwargs.get('aFe',0.0)
		radvel = kwargs.get('Vrad',0.0)
		rotvel = kwargs.get('Vrot',0.0)
		inst_R = kwargs.get('Inst_R',32000.0)
		logR = kwargs.get('logR',0.0)
		Dist = kwargs.get('Dist',100.0)
		Av = kwargs.get('Av',0.05)

		self.specpars = [Teff,logg,FeH,aFe,radvel,rotvel,inst_R]
		self.photpars = [Teff,logg,FeH,logR,Dist,Av]

		self.specmod = self.GM.genspec(self.specpars)
		self.sedmod  = self.GM.genphot(self.photpars)

	def writeout(self,*args,**kwargs):
		outfile_root = kwargs.get('output','test')
		overwrite = kwargs.get('overwrite',False)

		if overwrite:
			outfile = h5py.File('{}.h5'.format(outfile_root),'w')
			modindex = 0
		else:
			outfile = h5py.File('{}.h5'.format(outfile_root),'a')
			modindex = len(outfile.keys())

		outfile.create_dataset('mod_{0}/pars/specpars'.format(modindex),data=np.array(self.specpars))
		outfile.create_dataset('mod_{0}/pars/photpars'.format(modindex),data=np.array(self.photpars))
		outfile.create_dataset('mod_{0}/spec/wave'.format(modindex),data=np.array(self.specmod[0]))
		outfile.create_dataset('mod_{0}/spec/flux'.format(modindex),data=np.array(self.specmod[1]))

		for kk in self.sedmod.keys():
			outfile.create_dataset('mod_{0}/phot/{1}'.format(modindex,kk),data=self.sedmod[kk])

		outfile.close()


class SampleMIST(object):
	"""docstring for SampleMIST"""
	def __init__(self,):
		super(SampleMIST, self).__init__()
		MISTh5 = h5py.File('/Users/pcargile/Astro/MIST/MIST_v1.2/EEPTRACKS/MIST_1.2_EEPtrk_MSRGB.h5','r')
		self.MIST = {kk:MISTh5[kk] for kk in MISTh5['index']}
		self.gm = getmod()

	def pullsolar(self,):
		misttab = Table(np.array(self.MIST[b'0.00/0.00/0.00']))
		return misttab[misttab['initial_mass'] == 1.0]

	def mksolar(self,):
		soltab = self.pullsolar()
		for EEP in [300,450,600]:
			soltab_i = soltab[soltab['EEP'] == EEP]
			self.gm.gmod(Teff=10.0**soltab_i['log_Teff'][0],logg=soltab_i['log_g'][0],FeH=soltab_i['[Fe/H]'][0],aFe=0.0)
			self.gm.writeout(output='soltest')

if __name__ == '__main__':
	SM = SampleMIST()
	SM.mksolar()