#KJ66 Params: 
kj66_inputs = {
    'casing_od_inch': 4.33, #110mm
    'wall_thickness_mm': 0.5,
    'pressure_ratio': 2.2, #typical op. point (max ~2.4)
    'compressor_efficiency': 0.74, #cons est
    'mass_flow_air_kg_s': 0.23, #.22-.24
    'target_tit_k': 1123.0 #K
}
from V5_CombustionChamber_Design import *
kj66 = MicroJetCombustor(kj66_inputs)
results = kj66.run()
print_report(results, kj66_inputs)

#kj66_Actual