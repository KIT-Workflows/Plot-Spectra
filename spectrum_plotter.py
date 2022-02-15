import os, sys, glob, yaml

import numpy as np
from scipy.constants import c, h, hbar, e, m_e, N_A, epsilon_0
from scipy.stats import norm

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FuncFormatter

#from rdkit import Chem
#from rdkit.Chem import Draw
#from rdkit.Chem import MolFromSmiles as smi2mol

def Hart2eV(energy_H):
    return energy_H*m_e*(e**2/(4*epsilon_0*np.pi*hbar))**2/e    

def get_spectrum(transitions, x_range, smearing = 0.1):
    xvals = np.linspace(x_range[0],x_range[1],abs(int(100*(x_range[1]-x_range[0])+1)))
    yvals = np.zeros_like(xvals)

    for transition in transitions:
        if float(transition[1]) >= 0:
            yvals += norm.pdf((xvals-float(transition[0]))/smearing)*float(transition[1])

    return [xvals,yvals]

if __name__ == '__main__':
         
    with open('rendered_wano.yml') as infile:
        wano_file = yaml.full_load(infile)

    if wano_file['Plot IR spectra']:
        ir_files = glob.glob('*ir_result.yml')
        ir_files.sort()
        freq_range = [4000,400]

        for i in range(len(ir_files)):
            
            vibrations = []

            with open(ir_files[i]) as infile:
                ir_dict = yaml.full_load(infile)

            vib_list = ir_dict['vibrational frequencies']

            for vib in vib_list:
                if vib['frequency'] > 0:
                    vibrations.append([vib['frequency'],vib['intensity']])

            ir_spectrum = get_spectrum(vibrations, freq_range, smearing = 5)

            plt.gcf().clear
            fig, ax = plt.subplots()
            ax.set_xlabel('Wavenumber [cm^-1]')
            ax.set_ylabel('Intensity [a.u.]')
            ax.set_xlim(freq_range[0],freq_range[1])
    
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
          
            plt.plot(ir_spectrum[0], ir_spectrum[1], 'blue')

            if len(ir_files) == 1:
                plotname = 'IR_spectrum.png'

            else: 
                #plotname = '%s_IR_spectrum.png'%ir_dict['title']
                plotname = '%i_IR_spectrum.png'%i

            plt.savefig(plotname, dpi=300)

    
    absorptions, emissions = [], []
    abs_spectrum, em_spectrum = [], []

    nm_range = [300,700]
    ev_range = [(c*h*1E9/e)/nm for nm in nm_range]
    ev_range.sort()

    if wano_file['Plot UV/Vis absorption spectra']:
        abs_files = glob.glob('*abs_result.yml')
        abs_files.sort()

        if wano_file['Plot UV/Vis emission spectra']:
            em_files = glob.glob('*em_result.yml')
            em_files.sort()

        for i in range(len(abs_files)):

            with open(abs_files[i]) as infile: 
                abs_dict = yaml.full_load(infile)

            abs_list = abs_dict['electronic excitations']

            for abso in abs_list:
                absorptions.append([abso['exc energy'],abso['osc. strength']])

            for abso in absorptions:
                abso[0] = Hart2eV(abso[0])
    
            abs_spectrum = get_spectrum(absorptions, ev_range)

            if wano_file['Plot UV/Vis emission spectra']:
        
                with open(em_files[i]) as infile: 
                    em_dict = yaml.full_load(infile)

                em_list = em_dict['electronic excitations']

                for em in em_list:
                    emissions.append([em['exc energy'],em['osc. strength']])

                for em in emissions:
                    em[0] = Hart2eV(em[0])
    
                em_spectrum = get_spectrum(emissions, ev_range)

            if len(absorptions) > 0 or len(emissions) > 0:

                energy_ev = abs_spectrum[0]

                wavelength_nm = (c*h*1e9/e)/energy_ev

                plt.gcf().clear
                fig, ax = plt.subplots()
                ax.set_xlabel('Wavelength [nm]')
                ax.set_ylabel('Absorption [a.u.]')
                ax.set_xlim(nm_range[0],nm_range[1])
    
                ax.xaxis.set_minor_locator(AutoMinorLocator(2))
                ax.yaxis.set_minor_locator(AutoMinorLocator(2))
          
                top_ax=plt.twiny()
                top_ax.set_xlabel('Energy [eV]')
                top_ax.set_xlim(ax.get_xlim())
                top_ax.get_xaxis().set_major_formatter(FuncFormatter(lambda wavelength_nm, pos: "{:<04.3}".format((c*h*1E9/e)/wavelength_nm)))
                top_ax.xaxis.set_minor_locator(AutoMinorLocator(2))

                if len(abs_spectrum) > 1:
                    plt.plot(wavelength_nm, abs_spectrum[1],'blue')

                if len(em_spectrum) > 1:
                    plt.plot(wavelength_nm, em_spectrum[1],'red')

                plt.legend(['Absorption spectrum','Emission spectrum'],loc='upper right')

                #uv_dict = abs_dict

                if len(abs_files) == 1:
                    plotname = 'UV-Vis_spectrum.png'

                else: 
                    #plotname = '%s_UV-Vis_spectrum.png'%uv_dict['title']
                    plotname = '%i_UV-Vis_spectrum.png'%i

                plt.savefig(plotname, dpi=300)

    if wano_file['Delete yml files']: 
        for filename in glob.iglob('*ion_spectrum'):
            os.remove(filename)
