from dataclasses import dataclass
import numpy as np 
@dataclass
class genotype:
    x :int 
    y: int 
    x_bar: int
    y_bar: int

    def __post_init__(self):
        self.w = self.x + self.y 
        self.z = max(self.x_bar, self.y_bar)
        if self.w ==0:
            self.vaf = np.NaN
        else:
            self.vaf = self.z/self.w
    
    def __str__(self):
        return str(self.to_tuple())
    
    def __eq__(self, _value: object) -> bool:
        if isinstance(_value, genotype):
            return (self.x == _value.x) and (self.y == _value.y) \
                (self.x_bar == _value.x_bar) and (self.y_bar == _value.y_bar)
    
    def cna_eq(self, _value: object) -> bool:
        if isinstance(_value, genotype):
            return (self.x == _value.x) and (self.y == _value.y)
    
    def to_CNAgenotype(self):
        return CNAgenotype(self.x, self.y)
    
    def to_tuple(self):
        return (self.x, self.y, self.x_bar, self.y_bar)
@dataclass
class CNAgenotype:
    x: int 
    y: int 

    def __post_init__(self):
        self.total = self.x + self.y
 

    def __eq__(self, _value: object) -> bool:
        if isinstance(_value, CNAgenotype):
            return (self.x == _value.x) and (self.y == _value.y)
    
    def __lt__(self, _value: object) -> bool:
        if isinstance(_value, CNAgenotype):
            return self.total < _value.total
    
    def to_tuple(self):
        return (self.x, self.y)


# class genotype_list:
#     geno_list: list = geno_list




