#%%
from WGSPuzzler import WGSPuzzler
p=WGSPuzzler()
p.load_fasta("/Users/clavijob/sdg_py/cachn1/SC5314.fasta",[5])
p.make_mutated_copies(.01)
#p.make_mutated_copies(.0025)
p.pe_sequencing(error_ratio=0.02)
p.se_sequencing(coverage=25)
p.write_reference("puzzle_B-seed1_ref.fasta")
p.execute(prefix="puzzle_B-seed1",seed=1)
# %%
