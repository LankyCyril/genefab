from ._dataset import GLDS

class MicroarrayExperiment():
    """Implements wrapper of GLDS class that has 'Array Data Files'"""
    glds = None
    design_ref = None
    raw_data = None
    derived_data = None
    factors = None
 
    def __init__(self, glds):
        self.glds = glds
        self.factors = glds.factors()
        if glds.field_ids("Array Design REF"):
            self.design_ref = glds.property_table("Array Design REF")
        if glds.field_ids("Array Data File"):
            self.raw_data = glds.property_table("Array Data File")
        if glds.field_ids("Derived Array Data File"):
            self.derived_data = glds.property_table("Derived Array Data File")
        if (self.raw_data is None) and (self.derived_data is None):
            raise ValueError("No raw or derived data associated with factors")
