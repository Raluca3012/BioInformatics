# The melting temperature (Tm) is the temp at which 1/2 of a particular
# DNA seq. will dissociate and become a single strand of DNA. Primer length & seq.
# are of critical importance in the design of the params. over successful amplification
# (PCR). The melting temp. of a nucleic adic duplex increases both wit its length
# A simple formula used for the calculation of tm is 
# Tm = 4 * (G + C) + 2 * (A + T)
# The actual Tm is influenced by the concentration of Mg2+ , K+ , and cosolvents. An alternative formula is:
# Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
# where Na+ is the concentration of the solution and has a value of 0.001 
# Input 6-12 letters

import math 

def basic_tm(S):
    Tm = 4 * (S.count('G') + S.count('C')) + 2 * (S.count('A') + S.count('T'))
    return Tm

def advanced_tm(S):
    na = 0.001
    length = len(S)
    gc_count = S.count('G') + S.count('C')
    gc_perc = (gc_count / length) * 100
    log_na = math.log10(na)
    tm = 81.5 + 16.6 * log_na + 0.41 * gc_perc - (600 / length)
    return round(tm, 2)

S = 'TCCAGACGACTA'

print('The basic formula for Tm [Tm = 4 * (G + C) + 2 * (A + T)] has the following result:', basic_tm(S), '°C')

print('\nThe advanced formula for Tm [Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length] has the following result:', advanced_tm(S), '°C')
