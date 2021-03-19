from pymatgen.core import Lattice, Structure


def vary_lattice_dimensions(input_structure, varyamount=0.10, a=True, b=True, c=False):
    """ a tool to vary the lattice dimensions for a completed vasp run with the intent to check ISIF=3 and therefore
    not be so paranoid"""

    structure = Structure.from_file(input_pos)
    structure_new = structure.copy()
    a_var = 1
    b_var = 1
    c_var = 1
    if a:
        a_var += varyamount
    if b:
        b_var += varyamount
    if c:
        c_var = + varyamount
    structure_new.lattice = Lattice.from_parameters(a=a_var * structure.lattice.a,
                                                    b=b_var * structure.lattice.b,
                                                    c=c_var * structure.lattice.c,
                                                    alpha=structure.lattice.alpha,
                                                    beta=structure.lattice.beta,
                                                    gamma=structure.lattice.gamma)
    return structure_new
