# BabrahamLinkON
Analysis pipeline for VDJ-seq


conda install igblast



export BOWTIE2_INDEXES='/bi/scratch/Genomes/Mouse/GRCm38'
export BOWTIE2_REF='Mus_musculus.GRCm38'
#Home directory
export home='/bi/group/corcoran/Peter'
#Folder for all the log/output files
export log_folder=${home}/logs

#matplotlib backend for headless nodes
export MPLBACKEND=pdf

#specify tmp dir (needed for nodes as they don't have much memory)
export TMPDIR='/state/partition1'
