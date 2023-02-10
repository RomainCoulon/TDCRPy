import radioactivedecay as rd

Mo99_t0 = rd.Inventory({'Mo-99': 2.0}, 'Bq')
Mo99_t1 = Mo99_t0.decay(20.0, 'h')
Mo99_t1.activities('Bq')
