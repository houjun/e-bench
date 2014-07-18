ECoG-bench
=======

HDF5 read benchmark simulating I/O behavior of ECoG dataset analytics.

Group '/Descriptors'
Dataset 'ECoG_Fs'
(Sampling rate of ECoG data)
Dataset 'ELoc_CMF'
(vector of electrodes associated with area CMF of the brain)
Dataset 'ELoc_MTG'
(vector of electrodes associated with area MTG of the brain)
Dataset 'ELoc_SMG'
(vector of electrodes associated with area SMG of the brain)
Dataset 'ELoc_STG'
(vector of electrodes associated with area STG of the brain)
Dataset 'ELoc_parsO'
(vector of electrodes associated with area parsO of the brain)
Dataset 'ELoc_postCG'
(vector of electrodes associated with area postCG of the brain)
Dataset 'ELoc_preCG'
(vector of electrodes associated with area preCG of the brain)

1756 events, each with one trial type
Use 'Event_ECoGIndx' to find actual offsets in data, 'Event_EIndx' identifies the trial type of the event( and Event_ELbls is the name).

Dataset 'Event_ECoGIndx'
Size: 2x1756 (array onsets and offsets for the different ’trials’)

Dataset 'Event_EIndx'
(vector of indices for trial type)

Dataset 'Event_ELbls'
(string of trial names associated with indices in ‘Event_EIndx’
