from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='4uo0', model_segment=('FIRST:B','LAST:B'))
aln.append_model(mdl, align_codes='4uo0B', atom_files='4uo0.ent')
aln.append(file='2017-2018h3n2_HA_north', align_codes='ALL',alignment_format='FASTA')
aln.align2d()
aln.write(file='QNT-4uo0B.ali', alignment_format='PIR')
aln.write(file='QNT-4uo0B.pap', alignment_format='PAP')
