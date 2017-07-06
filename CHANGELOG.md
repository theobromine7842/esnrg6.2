2017-05-08 Kelly Tran

esnrg6.0.f
* changed 'TIP3' to 'TP4E'
* changed 'NA' to 'POT'
* changed 'read' for:
	* calc_type='ESE'
	* topfile='top_all36_prot.rtf'
	* paramfile='par_all36_prot.prm'
	* cut_type='none'
	* num_dcd=1
	* PBC_type='cube'
	* rmin=1
	* rmax=15
	* div=1
	* SOLV_quest='N'
	* strtframe=icntrl(1)
	* endframe=icntrl(1)+icntrl(3)
	* skipframe=icntrl(2)
* changed:
	* totframe = icntrl(3)
	* write(6,'(A,I8)') 'Number of frames = ',totframe
	* write(6,'(A,I9)') 'Initial step = ',icntrl(1)
	* write(6,'(A,I6)') 'Write frequency = ',icntrl(2)

	
* write data files to dat/ directory
* added "numframe = 0" to label 206


assigngrp.f
* changed 'TIP' to 'TP4'
* added "if (atid.eq.'POT') gtype(numgrp) = 6" and "if (atid.eq.'CLA') gtype(numgrp) = 6"


scp knt8@medusa.georgetown.edu:/home/knt8/3-cafd/02-anal-ions/parabolas/esnrg6.0/esnrg6.0.f .


Number of frames in this file =      100
Number of previous steps =    101000
Frequency (steps) for saving frames =      0
