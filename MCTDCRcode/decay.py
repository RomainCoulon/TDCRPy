import radioactivedecay as rd

Ac225_t0 = rd.Inventory({'Ac-225': 1}, 'Bq')
Ac225_t1 = Ac225_t0.decay(2, 'd')
Ac225_t1.activities('Bq')
print(Ac225_t1)
