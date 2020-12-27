import os
import numpy as np
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
import matplotlib.pyplot as plt
    
def precise_align(env,tmp,targ,tmp_pdb,model_segment=None):
    """use `align2d` function to align target and template using dynamic programming"""
    tmp_pdb = os.path.relpath(tmp_pdb) #get relative path to avoid colon errors
    aln = alignment(env)
    mdl = model(env, file=tmp, model_segment=model_segment)
    tmp = tmp+''.join(set((model_segment[0][-1]+model_segment[1][-1])))
    aln.append_model(mdl, align_codes=tmp, atom_files=tmp_pdb)
    aln.append(file=targ+'.ali', align_codes=targ)
    aln.align2d()
    
    aln.write(file=targ+'-'+tmp+'.ali', alignment_format='PIR')
    aln.write(file=targ+'-'+tmp+'.pap', alignment_format='PAP')
    
def ali_to_faa(file):
    """make a `.faa` file from a `.ali` file"""
    sequence_data = open(file).readlines(); sequence_data.pop(1) #unwinds into many lines and removes the second line
    sequence_data = np.sum(np.array(sequence_data,dtype='object')) #stitches the string back
    faa_file = open(file[:-4]+'.faa', "w")
    faa_file.write(sequence_data); faa_file.close()
    
def model_maker(env,tmp,targ):
    """Generate a structure model using satisfaction of spatial restraints"""
    a = automodel(env, alnfile=targ+'-'+tmp+'.ali',
                  knowns=tmp, sequence=targ,
                  assess_methods=(assess.DOPE,
                                  assess.GA341))
    a.starting_model = 1
    a.ending_model = 5
    a.make()
    
def dope_evaluate(env,file_name,profile_name,model_segment=None):
    """Evaluate the `.pdb` file and generate a DOPE profile"""
    env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
    env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

    # read model file
    mdl = complete_pdb(env, file_name, model_segment=model_segment)

    # Assess with DOPE:
    s = selection(mdl)   # all atom selection
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=profile_name,
                  normalize_profile=True, smoothing_window=15)

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in reversed(list(enumerate(seq.residues))):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals

def dope_plot(env,tmp,targ,save=False):
    """Plot the dope profile of the model and the template protein"""
    a = alignment(env, file=targ+'-'+tmp+'.ali')
    template = get_profile(tmp+'.profile', a[tmp])
    model = get_profile(targ+'.profile', a[targ])

    # Plot the template and model profiles in the same plot for comparison:
    plt.figure(figsize=(7,4))
    plt.xlabel('Alignment Position')
    plt.ylabel('DOPE per-residue score')
    plt.plot(model, linewidth=2, label='Model')
    plt.plot(template, linewidth=2, label='Template')
    plt.legend()
    if save:
        plt.savefig('dope_profile.png', dpi=300)
    plt.show()