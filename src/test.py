#%%
from WGSPuzzler import WGSPuzzler
p=WGSPuzzler()
p.load_fasta("/Users/clavijob/sdg_py/cachn1/SC5314.fasta",[5])

# p.make_mutated_copies(.01)
# p.make_mutated_copies(.0025)
# p.pe_sequencing(error_ratio=0.02)
# p.se_sequencing(coverage=25)
# p.write_reference("puzzle_A-seed1_ref.fasta")
# p.execute(prefix="puzzle_A-seed1",seed=1)

# p.make_mutated_copies(.01)
# p.pe_sequencing(error_ratio=0.02)
# p.se_sequencing(coverage=25)
# p.write_reference("puzzle_B-seed1_ref.fasta")
# p.execute(prefix="puzzle_B-seed1",seed=1)

# p.make_mutated_copies(.01)
# p.make_mutated_copies(.01)
# p.pe_sequencing(error_ratio=0.02)
# p.se_sequencing(coverage=25)
# p.write_reference("puzzle_C-seed1_ref.fasta")
# p.execute(prefix="puzzle_C-seed1",seed=1)

#p.make_mutated_copies(.01)
#p.make_mutated_copies(.02)
#p.pe_sequencing(error_ratio=0.02)
#p.se_sequencing(coverage=25)
#p.write_reference("puzzle_D-seed1_ref.fasta")
#p.execute(prefix="puzzle_D-seed1",seed=1)

p.make_mutated_copies(.01,[0])
p.make_mutated_copies(.01,[0])
p.select_molecules([1,2])
p.make_mutated_copies(.02)
p.pe_sequencing(error_ratio=0.02)
p.se_sequencing(coverage=25)
p.write_reference("puzzle_E-seed1_ref.fasta")
p.execute(prefix="puzzle_E-seed1",seed=1)
# %%
