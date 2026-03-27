import math

class MicroJetCombustor:
    def __init__(self, inputs):
        self.inputs = inputs
        self.res = {}
        self.R = 287.05     # J/(kg·K) — specific gas constant for air
        self.GAMMA = 1.4    # ratio of specific heats for air

        self.DESIGN_PARAMS = {
            # --- Annulus Velocities ---
            'target_annulus_vel':       35.0,  # m/s — OUTER annulus feed velocity
            #   Higher velocity needed here: outer path carries 6-12 vaporizer scoops that
            #   require momentum to fill cleanly. Too slow → scoop starvation → rich primary.
            'target_inner_annulus_vel': None,  # m/s — INNER annulus feed velocity target.
            #   None (default) = auto-match outer annulus velocity.
            #   This is the physically correct general default: both annuli carry air at
            #   the same bulk velocity, giving equal static pressure recovery on both sides
            #   of the liner and symmetric hydraulic resistance around the combustion zone.
            #
            #   The previous hard-coded 20 m/s was an empirical quirk of the KJ66 geometry
            #   (small engine, thick inner annulus relative to combustion gap) that produced
            #   correct KJ66 dimensions by construction — but poisoned scaling for any other
            #   engine size. At 0.487 kg/s (6-inch engine), 20 m/s inflates the inner annulus
            #   to 27.5 mm, crushing the combustion gap down to 11.4 mm despite being in a
            #   152 mm casing. The fix: default to v_outer, override in presets that need it.
            #
            #   Override options:
            #     None        → match v_outer (correct general default)
            #     20.0        → KJ66 empirical calibration (set in kj66_inputs preset)
            #     any float   → explicit target for custom designs
            'target_tube_liq_vel': 3.0,     # m/s  — DEPRECATED: superseded by target_vap_mix_vel
            # Vaporizer bore target: gas-phase air+fuel mixture velocity through the bore.
            # At ~85 m/s and 600K, 1.5 g/s per tube gives ID ≈ 4.2 mm — matching the KJ66 empirical
            # 5mm OD / 4mm ID tube stock. Correct physics: the bore carries the gas-dominant mixture,
            # not pure liquid fuel. Range: 70–100 m/s.
            'target_vap_mix_vel': 85.0,     # m/s  — gas-phase mixture velocity in vaporizer bore
            'vap_pitch_mm': 50.0,           # mm   — circumferential pitch between vaporizer scoops
            # Physical basis: KJ66 uses 6 tubes at 60° spacing around the outer liner ID
            # (π × 95.8 mm ÷ 6 ≈ 50.2 mm). Previously 35 mm was used with mean combustion
            # diameter as the reference — that combination gave n_raw ≈ 7.1 → snapped to 8.
            # With the dome-rim (outer_liner_id) reference and 50 mm pitch: n_raw ≈ 6.0 → 6. ✓
            # --- Combustor Sizing ---
            'target_pressure_drop': 0.04,   # 4%   — dP across liner (industry standard for micro-GT)
            'discharge_coeff_hole': 0.60,   # Cd   — sharp-edged holes (Lefebvre recommendation)
            # --- Zone Length Fractions of Total Liner Length ---
            'frac_primary_zone': 0.30,      # 30%  — primary/rich zone (flame stabilization)
            'frac_secondary_zone': 0.30,    # 30%  — secondary zone (soot burnout, CO oxidation)
            'frac_dilution_zone': 0.40,     # 40%  — dilution zone (temperature profile shaping)
            # --- Feed Annulus Air Split ---
            'f_outer_feed': 0.60,
            # --- Combustor Length Minimums ---
            # tau_min: minimum residence time floor for chamber length sizing (seconds).
            #   0.002 (2ms) reproduces KJ66 empirical L/D≈0.70 (~55mm for KJ66 geometry).
            #   Modern designs: use 0.003–0.004 for better combustion completeness.
            'tau_min_s': 0.002,
            # L_D_min: minimum chamber length as multiple of mean combustion diameter.
            #   0.70 reproduces KJ66 geometry. Modern annular designs: 0.90–1.2.
            'L_D_min': 0.70,
            # --- Hole Split Between Inner and Outer Liner ---
            'outer_liner_hole_fraction': None,  # e.g. 0.60; None = auto from circumference ratio
            # --- Air Budget Fractions ---
            'film_cooling_fraction': 0.12,  # 12%  — reserved for liner film cooling / effusion
            # --- Vaporizer Air Ingestion ---
            'primary_air_vaporizer_fraction': 0.12,
            # --- Vaporizer Tube Count Constraint ---
            'allowed_vaporizer_counts': [6, 8, 12],
            # --- Vaporizer Scoop Inlet ---
            'vap_scoop_target_vel_m_s': 22.0,  # m/s — low velocity for stable entrainment
            'vap_scoop_cd':             0.72,  # Cd — flared bell-mouth or open tube
            # --- Explicit Film-Cooling Pattern ---
            'film_row_pitch_mm':    38.0,  # mm  — axial spacing between cooling rows
            'film_hole_dia_mm':      1.20,  # mm  — discrete effusion hole diameter
            'film_rows_primary':        1,  # rows in primary zone
            'film_rows_secondary':      1,  # rows in secondary zone
            'film_rows_dilution':       2,  # rows in dilution zone (longest, most exposed)
            # --- Pressure-Balance Solver Loss Coefficients ---
            'K_turn':     1.0,   # 180° dome turn loss, outer path
            'K_entrance': 0.5,   # sharp-edged inner annulus entry (Idelchik §5-1)
            # --- Fuel System Target Pressure Drop ---
            'fuel_dP_target_bar': 0.60,   # bar — target total fuel system ΔP at full thrust
            # --- Secondary Zone Hole Count Multiplier ---
            # 2 holes per vaporizer per liner (12 holes total at n_vap=6).
            # At 1×/vap: d_sec_out ≈ 8.7 mm → Sturgess penetration 1.44×H (inner-liner impingement).
            # At 2×/vap: d_sec_out ≈ 6.1 mm → penetration 1.01×H (acceptable with dual-wall coverage).
            # Reference: empirical KJ66 secondary hole diameters 4–6 mm outer [S5].
            'sec_holes_per_vap': 2,
            # --- Dilution Zone Hole Count Multiplier ---
            # Decoupled from vaporizer count per Reviewer 1 recommendation.
            # KJ66 empirical: 2×/vap (12 holes/liner) — reproduces historical geometry.
            # Better mixing: 4×/vap (24 holes/liner) — reduces penetration depth ~30%.
            # Custom engines default to 4; set to 2 in the KJ66 preset for historical fidelity.
            'dil_holes_per_vap': 4,
        }

        # --- Fuel Properties (Jet-A / Kerosene) ---
        self.FUEL = {
            "LHV":       43.0e6,    # J/kg  — lower heating value
            "STOICH_AFR": 14.7,     # kg_air/kg_fuel — stoichiometric air-fuel ratio
            "RHO_LIQ":   800.0,     # kg/m³ — liquid density at ambient
            "T_BOIL":    450.0,     # K     — approximate vaporization temperature
        }

        self.MATERIALS = {
            '304SS':  {'alpha': 17.2e-6, 'name': '304 Stainless Steel'},
            '316SS':  {'alpha': 16.0e-6, 'name': '316 Stainless Steel'},
            'IN625':  {'alpha': 13.0e-6, 'name': 'Inconel 625'},
            'IN718':  {'alpha': 13.0e-6, 'name': 'Inconel 718'},
        }

    def _get_cp(self, T_kelvin):
        """
        Variable specific heat for air/combustion products. Quadratic fit vs T.

        Calibrated against NACA SP-273 / Gordon-McBride equilibrium data for
        lean-to-stoichiometric Jet-A combustion products at 1-4 atm:
          300K→1005, 600K→1050, 900K→1090, 1200K→1130,
          1500K→1175, 1800K→1250, 2100K→1290, 2200K→1310  J/kgK

        Quadratic: cp = 1006.9 + 0.12415*(T-300) + 1.919e-5*(T-300)^2
        Max fit error: ±14 J/kgK at 1800K (~1.1%) — acceptable for 1D sizing.

        Previous model (1005 + 0.10T, capped at 1150) under-predicted cp by
        100–160 J/kgK above 1800K and was non-physical above 2000K.
        The quadratic is monotonically increasing and unbounded (physically correct
        for lean products; for richer zones the model is conservative).
        """
        T = max(300.0, min(float(T_kelvin), 2200.0))
        return 1006.942 + 0.124151 * (T - 300.0) + 1.9189e-5 * (T - 300.0) ** 2

    def thermodynamics(self):
        P_amb = 101325
        T_amb = 288.15
        PR    = self.inputs['pressure_ratio']
        eta_c = self.inputs['compressor_efficiency']
        P2     = P_amb * PR
        T2_iso = T_amb * (PR ** ((self.GAMMA - 1) / self.GAMMA))
        T2     = T_amb + (T2_iso - T_amb) / eta_c
        rho2   = P2 / (self.R * T2)
        self.res['P2_Pa']  = P2
        self.res['T2_K']   = T2
        self.res['rho2']   = rho2

    def mass_flow_and_fuel(self):
        """Determines air mass flow and required fuel flow for target TIT."""
        m_air = self.inputs.get('mass_flow_air_kg_s')
        if m_air is None:
            ref_od   = 0.1524
            ref_flow = 0.45
            od_m     = self.inputs['casing_od_inch'] * 0.0254
            m_air    = ref_flow * (od_m ** 2 / ref_od ** 2)

        target_tit = self.inputs['target_tit_k']
        T2         = self.res['T2_K']
        T_avg      = (target_tit + T2) / 2
        cp_avg     = self._get_cp(T_avg)
        comb_eff   = 0.96

        # heat_loss_factor: fraction of released combustion energy that reaches the
        # working gas vs escaping through the casing walls.
        #
        # CORRECT THERMODYNAMIC DIRECTION (Rev 15 fix):
        #   m_fuel * LHV * eta_comb * eta_loss = m_air * cp * (TIT - T2)
        #   → m_fuel = m_air * cp * (TIT - T2) / (LHV * eta_comb * eta_loss)
        #   With eta_loss < 1: MORE fuel is required (some heat escapes, so the
        #   combustor must burn more to still reach TIT). This is physically correct.
        #
        # Rev 13/14 incorrectly applied the factor as a multiplier on delta_T
        # (effective_dT = dT * hlf), which made less heat-loss require LESS fuel —
        # thermodynamically backwards. That version calibrated by accident: the
        # 0.88 multiplier happened to produce ~3.8 g/s by compressing the energy
        # demand rather than increasing the denominator. The result was plausible
        # numerically but physically wrong.
        #
        # Correct calibration: adiabatic (eta_loss=1.0) gives mdot_fuel=4.38 g/s
        # at KJ66 full-thrust TIT=1123K. Schreckling's 3.8 g/s is part-throttle
        # (~1040K TIT equivalent) — not a valid calibration target for full-thrust
        # design-point sizing. The adiabatic value is the correct design-point fuel
        # flow when sizing hole areas, zone lengths, and vaporizer flow rates.
        #
        # Override range: 0.80–1.00 (1.00 = default adiabatic; 0.85 = upper bound
        # for small engines with uninsulated casings, per Rodgers 1991)
        heat_loss_factor = self.inputs.get('heat_loss_factor', 1.0)
        energy_req = m_air * cp_avg * (target_tit - T2)
        m_fuel     = energy_req / (self.FUEL['LHV'] * comb_eff * heat_loss_factor)

        self.res['mdot_air']          = m_air
        self.res['mdot_fuel']         = m_fuel
        self.res['overall_AFR']       = m_air / m_fuel
        self.res['overall_phi']       = self.FUEL['STOICH_AFR'] / self.res['overall_AFR']
        self.res['cp_used']           = cp_avg
        self.res['heat_loss_factor']  = heat_loss_factor

    def zonal_analysis(self):
        m_fuel  = self.res['mdot_fuel']
        stoich  = self.FUEL['STOICH_AFR']
        m_total = self.res['mdot_air']

        f_film  = self.DESIGN_PARAMS['film_cooling_fraction']
        m_film  = m_total * f_film
        m_comb  = m_total - m_film

        phi_pri       = 1.6
        m_air_pri_total = (m_fuel * stoich) / phi_pri

        if m_air_pri_total > m_comb:
            raise ValueError(
                f"Primary zone requires {m_air_pri_total*1000:.1f} g/s air but only "
                f"{m_comb*1000:.1f} g/s is available after film cooling reservation. "
                f"Reduce film_cooling_fraction, increase mass_flow_air, or lower target_tit_k."
            )

        f_vap         = self.DESIGN_PARAMS['primary_air_vaporizer_fraction']
        m_air_vap     = m_air_pri_total * f_vap
        m_air_pri_liner = m_air_pri_total * (1.0 - f_vap)
        m_air_pri = m_air_pri_total

        # phi_sec_target = 0.60: standard Lefebvre §4.3 value for annular combustors.
        # Rev13 raised this to 0.70 to rescue dilution air starved by excess fuel flow.
        # With heat_loss_factor fix reducing mdot_fuel to ~3.8 g/s, dilution is no longer
        # starved and phi_sec=0.60 is correct. At phi=0.60, CO oxidation is complete and
        # secondary zone exits ~1500-1600 K — correct entry temperature for dilution quench.
        phi_sec_target = 0.60
        m_cumulative   = (m_fuel * stoich) / phi_sec_target
        m_air_sec      = m_cumulative - m_air_pri

        m_air_dil = m_comb - m_air_pri - m_air_sec
        if m_air_dil < 0:
            m_air_dil = 0
            m_air_sec = m_comb - m_air_pri

        self.res['m_film_cooling']      = m_film
        self.res['split_primary']       = m_air_pri
        self.res['split_primary_liner'] = m_air_pri_liner
        self.res['split_primary_vap']   = m_air_vap
        self.res['split_secondary']     = m_air_sec
        self.res['split_dilution']      = m_air_dil
        self.res['m_combustion_air']    = m_comb

        # Actual primary zone equivalence ratio (using ALL primary air — liner + vaporizer)
        phi_primary_actual = m_fuel / (m_air_pri / stoich) if m_air_pri > 0 else 0.0
        self.res['phi_primary_actual']  = phi_primary_actual

        cp_flame        = self._get_cp(1900)
        m_fuel_burned   = m_air_pri / stoich
        Q_zone1         = m_fuel_burned * self.FUEL['LHV']
        m_zone1_total   = m_air_pri + m_fuel
        T_pri_raw       = self.res['T2_K'] + Q_zone1 / (m_zone1_total * cp_flame)
        T_pri_clamped   = min(T_pri_raw, 1950.0)
        # Cap at 1950 K: at φ_pri=1.6 (rich), incomplete combustion and radiation
        # losses limit the effective bulk primary zone temperature to ~1850–1950 K.
        # The raw adiabatic calculation produces higher values (φ=1.6 flame temperature
        # for Jet-A is ~2100 K) because it ignores the unburned fuel fraction and
        # radiation losses. 1950 K is the appropriate screening cap for a rich primary
        # zone; 2050 K was the lean/stoichiometric adiabatic limit (wrong zone type).
        # Impact: dome heat flux and J_primary calculations use T_primary — a lower cap
        # reduces h_conv and rho_hot, giving slightly more conservative J values.
        self.res['T_primary_zone_est']     = T_pri_clamped
        self.res['T_primary_zone_raw']     = T_pri_raw
        self.res['T_primary_zone_clamped'] = T_pri_raw > 1950.0   # matches the cap on the line above

    def mechanical_geometry(self):
        wall_m     = self.inputs['wall_thickness_mm'] / 1000.0
        od_casing  = self.inputs['casing_od_inch'] * 0.0254
        id_casing  = od_casing - 2 * wall_m
        od_tunnel  = self.inputs['shaft_tunnel_od_inch'] * 0.0254
        id_tunnel  = od_tunnel - 2 * wall_m

        rho2   = self.res['rho2']
        m_air  = self.res['mdot_air']
        v_outer    = self.DESIGN_PARAMS['target_annulus_vel']
        # Resolve inner annulus velocity — three-tier priority:
        #   1. inputs dict   (per-run override, e.g. kj66_inputs['target_inner_annulus_vel'] = 20.0)
        #   2. DESIGN_PARAMS (class default, now None — skip)
        #   3. Fall back to v_outer  (symmetric velocity = correct general default)
        # This preserves the KJ66 empirical 20 m/s when explicitly set in the inputs dict,
        # while all other engine scales default to equal-velocity annuli (correct physics).
        _dp_inner = self.DESIGN_PARAMS.get('target_inner_annulus_vel')
        v_inner   = (self.inputs.get('target_inner_annulus_vel')
                     or (_dp_inner if _dp_inner is not None else None)
                     or v_outer)
        f_outer_feed = self.DESIGN_PARAMS['f_outer_feed']

        A_outer_feed = (m_air * f_outer_feed) / (rho2 * v_outer)
        A_inner_feed = (m_air * (1.0 - f_outer_feed)) / (rho2 * v_inner)

        r_casing_id    = id_casing / 2
        A_casing_id    = math.pi * r_casing_id ** 2
        if A_outer_feed >= A_casing_id:
            raise ValueError(
                f"Outer feed annulus area ({A_outer_feed*1e6:.0f} mm²) exceeds casing bore "
                f"area ({A_casing_id*1e6:.0f} mm²). Mass flow or casing OD is infeasible: "
                f"reduce mass_flow_air, increase casing_od_inch, or increase target_annulus_vel."
            )
        r_outer_liner_od = math.sqrt(r_casing_id ** 2 - A_outer_feed / math.pi)
        outer_liner_od = r_outer_liner_od * 2

        r_tunnel_od    = od_tunnel / 2
        r_inner_liner_id = math.sqrt(r_tunnel_od ** 2 + A_inner_feed / math.pi)
        inner_liner_id = r_inner_liner_id * 2

        # Guard: inner liner ID must clear the shaft tunnel OD with a minimum annulus gap.
        # The previous check (r_tunnel_od**2 + A_inner_feed/π < 0) was mathematically
        # impossible — both terms are non-negative. The real failure mode is the computed
        # inner_liner_id falling at or below od_tunnel (zero or negative annulus gap).
        inner_liner_id_check = math.sqrt(r_tunnel_od ** 2 + A_inner_feed / math.pi) * 2
        if A_inner_feed <= 0 or inner_liner_id_check <= od_tunnel * 1.05:
            raise ValueError(
                f"Inner feed annulus is infeasible: computed inner liner ID "
                f"({inner_liner_id_check*1000:.1f} mm) is at or below shaft tunnel OD "
                f"({od_tunnel*1000:.1f} mm). "
                f"Reduce inner annulus velocity target, increase casing OD, or reduce shaft tunnel OD."
            )

        outer_liner_id = outer_liner_od - 2 * wall_m
        inner_liner_od = inner_liner_id + 2 * wall_m

        A_comb = math.pi * ((outer_liner_id / 2) ** 2 - (inner_liner_od / 2) ** 2)
        D_mean = (outer_liner_id + inner_liner_od) / 2
        combustion_gap = (outer_liner_id - inner_liner_od) / 2

        H_liner_span = (outer_liner_od - inner_liner_id) / 2
        gap_min = max(0.008, H_liner_span * 0.60)
        if combustion_gap <= gap_min:
            raise ValueError(
                f"Geometry infeasible: combustion annulus gap = {combustion_gap*1000:.1f} mm "
                f"(minimum {gap_min*1000:.1f} mm = 60% of liner span "
                f"{H_liner_span*1000:.1f} mm). "
                f"Reduce shaft tunnel OD, increase casing OD, or reduce wall thickness."
            )

        tau_min   = self.inputs.get('tau_min_s',
                        self.DESIGN_PARAMS.get('tau_min_s', 0.002))
        L_D_min   = self.inputs.get('L_D_min',
                        self.DESIGN_PARAMS.get('L_D_min',   0.70))
        m_comb    = self.res.get('m_combustion_air', m_air)
        # Use bulk combustion temperature for density — the gas in the chamber is hot.
        # Cold inlet density (rho2) understates V_ref by ~2×, making L_min_tau ~2× too short.
        # T_bulk = mean of inlet and exit temperatures — consistent with tau_comb calculation.
        T_bulk   = (self.res['T2_K'] + self.inputs['target_tit_k']) / 2.0
        rho_bulk = self.res['P2_Pa'] / (self.R * T_bulk)
        V_ref     = m_comb / (rho_bulk * A_comb)
        L_min_tau = V_ref * tau_min
        L_min_geom = D_mean * L_D_min
        length = max(L_min_tau, L_min_geom)

        m_air_vap = self.res.get('split_primary_vap', 0.0)
        A_outer_feed_actual = math.pi * (r_casing_id ** 2 - r_outer_liner_od ** 2)
        A_inner_feed_actual = math.pi * (r_inner_liner_id ** 2 - r_tunnel_od ** 2)
        m_outer_entry = m_air * f_outer_feed
        m_inner_entry = m_air * (1.0 - f_outer_feed)
        v_outer_actual = m_outer_entry / (rho2 * A_outer_feed_actual) if A_outer_feed_actual > 0 else 0
        v_inner_actual = m_inner_entry / (rho2 * A_inner_feed_actual) if A_inner_feed_actual > 0 else 0
        m_outer_post_vap = m_outer_entry - m_air_vap
        v_outer_post_vap = m_outer_post_vap / (rho2 * A_outer_feed_actual) if A_outer_feed_actual > 0 else 0

        self.res['f_outer_feed_used']     = f_outer_feed
        self.res['casing_od_mm']          = od_casing   * 1000
        self.res['casing_id_mm']          = id_casing   * 1000
        self.res['outer_liner_od_mm']     = outer_liner_od * 1000
        self.res['outer_liner_id_mm']     = outer_liner_id * 1000
        self.res['inner_liner_od_mm']     = inner_liner_od * 1000
        self.res['inner_liner_id_mm']     = inner_liner_id * 1000
        self.res['shaft_tunnel_od_mm']    = od_tunnel   * 1000
        self.res['shaft_tunnel_id_mm']    = id_tunnel   * 1000
        self.res['combustion_gap_mm']     = combustion_gap * 1000
        self.res['combustion_annulus_A']  = A_comb
        self.res['D_mean_comb_mm']        = D_mean     * 1000
        self.res['chamber_length_mm']     = length     * 1000
        self.res['outer_annulus_gap_mm']  = (id_casing - outer_liner_od) / 2 * 1000
        self.res['inner_annulus_gap_mm']  = (inner_liner_id - od_tunnel) / 2 * 1000
        self.res['v_outer_annulus']       = v_outer_actual
        self.res['v_inner_annulus']       = v_inner_actual
        self.res['v_outer_post_vap']      = v_outer_post_vap
        self.res['outer_liner_circumf_m'] = math.pi * outer_liner_id
        self.res['inner_liner_circumf_m'] = math.pi * inner_liner_od

        mat_key  = self.inputs.get('liner_material', '304SS').upper()
        if mat_key not in self.MATERIALS:
            raise ValueError(f"Unknown liner_material '{mat_key}'. Valid options: {list(self.MATERIALS.keys())}")
        mat      = self.MATERIALS[mat_key]
        alpha    = mat['alpha']
        delta_T  = self.inputs['target_tit_k'] - 293.0

        dR_outer = alpha * (outer_liner_od / 2) * delta_T * 1000
        dR_inner = alpha * (inner_liner_od / 2) * delta_T * 1000
        dR_outer_clamped = max(0.10, min(1.00, dR_outer))   # ceiling 0.30→1.00 (KJ66 outer: 0.691mm)
        dR_inner_clamped = max(0.10, min(1.00, dR_inner))   # ceiling 0.30→1.00 (KJ66 inner: 0.448mm)
        outer_liner_od_cold = outer_liner_od * 1000 - dR_outer_clamped
        inner_liner_od_cold = inner_liner_od * 1000 - dR_inner_clamped

        self.res['liner_material']              = mat_key
        self.res['liner_material_name']         = mat['name']
        self.res['liner_alpha_per_K']           = alpha
        self.res['outer_liner_od_hot_mm']       = round(outer_liner_od * 1000, 3)
        self.res['outer_liner_od_cold_mm']      = round(outer_liner_od_cold, 3)
        self.res['outer_liner_thermal_offset_mm'] = round(dR_outer_clamped, 3)
        self.res['outer_liner_thermal_offset_unclamped_mm'] = round(dR_outer * 1000, 4) / 1000  # for flag
        self.res['inner_liner_od_hot_mm']       = round(inner_liner_od * 1000, 3)
        self.res['inner_liner_od_cold_mm']      = round(inner_liner_od_cold, 3)
        self.res['inner_liner_thermal_offset_mm'] = round(dR_inner_clamped, 3)
        self.res['inner_liner_thermal_offset_unclamped_mm'] = round(dR_inner * 1000, 4) / 1000
        self.res['thermal_offset_warn_outer'] = dR_outer > 0.50   # flag: offset capped, significant growth
        self.res['thermal_offset_warn_inner'] = dR_inner > 0.50

        self.res['L_primary_mm']   = length * self.DESIGN_PARAMS['frac_primary_zone'] * 1000
        self.res['L_secondary_mm'] = length * self.DESIGN_PARAMS['frac_secondary_zone'] * 1000
        self.res['L_dilution_mm']  = length * self.DESIGN_PARAMS['frac_dilution_zone'] * 1000

    def vaporizer_tubes(self):
        dome_circumf_m = math.pi * self.res['outer_liner_id_mm'] / 1000
        pitch_m        = self.DESIGN_PARAMS['vap_pitch_mm'] / 1000
        n_raw          = dome_circumf_m / pitch_m
        mean_circumf_m = math.pi * self.res['D_mean_comb_mm'] / 1000

        allowed = self.DESIGN_PARAMS['allowed_vaporizer_counts']
        # Snap logic — two-stage:
        #   Stage 1 (near-integer): if n_raw is within 0.15 of an allowed count,
        #     use nearest-neighbour snap. This handles the common case where the
        #     pitch produces n_raw = 6.02 and we correctly stay at 6, not jump to 8.
        #   Stage 2 (between counts): use ceiling snap — smallest allowed n ≥ n_raw.
        #     Rationale: n_raw=6.5 means 6 tubes are spaced wider than the design
        #     pitch, reducing scoop coverage uniformity; 8 is the conservative choice.
        # The previous math.ceil(n_raw - 1e-9) logic converted 6.02→7→snap to 8. Wrong.
        snap_tol = 0.15
        near = [n for n in sorted(allowed) if abs(n - n_raw) <= snap_tol]
        if near:
            n_tubes = min(near, key=lambda n: abs(n - n_raw))   # nearest in tolerance band
        else:
            candidates = [n for n in sorted(allowed) if n >= n_raw]
            n_tubes = candidates[0] if candidates else max(allowed)  # strict ceiling

        m_air_vap       = self.res['split_primary_vap']
        m_fuel_total    = self.res['mdot_fuel']
        m_total_vap     = m_fuel_total + m_air_vap
        m_mix_per_tube  = m_total_vap / n_tubes
        m_fuel_per_tube = m_fuel_total / n_tubes

        # ── Tube bore from gas-phase mixture velocity ────────────────────────
        # The bore carries the air+fuel gas mixture at ~600 K.
        # Sizing from the target gas mixture velocity gives the correct bore
        # from physics, not from a minimum-floor band-aid.
        # At 85 m/s: KJ66 → 4.22 mm ID (matches empirical 5mm OD / 4mm ID stock).
        T_vap_mix    = 600.0   # K — conservative mixture temperature in bore
        rho_mix_bore = self.res['P2_Pa'] / (self.R * T_vap_mix)
        target_v_mix = self.DESIGN_PARAMS.get('target_vap_mix_vel', 85.0)
        A_bore       = m_mix_per_tube / (rho_mix_bore * target_v_mix)
        id_tube      = math.sqrt(4 * A_bore / math.pi)
        id_tube      = max(id_tube, 0.003)   # 3 mm floor: practical machining minimum for SS tubing
        od_tube      = id_tube + 0.0012      # ~1.2 mm wall

        # ── Tube friction ΔP — gas-phase Darcy-Weisbach ──────────────────────
        # Use gas mixture properties (rho at 600K, mu_air at 600K), not liquid fuel.
        mu_mix_bore = 3.0e-5   # Pa·s — air at 600 K (mixture approximation)
        Re_tube     = rho_mix_bore * target_v_mix * id_tube / mu_mix_bore
        if Re_tube < 2300:
            f_darcy = 64.0 / Re_tube
        else:
            f_darcy = 0.316 * Re_tube ** (-0.25)   # Blasius

        L_tube_eff = (self.res.get('chamber_length_mm', 120) / 1000) * 1.5 + 2 * od_tube
        dP_tube_Pa = f_darcy * (L_tube_eff / id_tube) * (rho_mix_bore * target_v_mix ** 2 / 2)

        T_exit_gas    = 950.0
        rho_gas_exit  = self.res['P2_Pa'] / (self.R * T_exit_gas)
        Cd_crimp        = 0.61
        dP_target_Pa    = self.DESIGN_PARAMS['fuel_dP_target_bar'] * 1e5
        dP_crimp_req_Pa = dP_target_Pa - dP_tube_Pa

        if dP_crimp_req_Pa > 0:
            v_crimp    = math.sqrt(2 * Cd_crimp**2 * dP_crimp_req_Pa / rho_gas_exit)
            A_crimp    = m_mix_per_tube / (rho_gas_exit * v_crimp)
            d_crimp_mm = math.sqrt(4 * A_crimp / math.pi) * 1000
            d_crimp_mm = max(d_crimp_mm, 0.8)
            A_crimp_actual = math.pi * (d_crimp_mm / 2000) ** 2
            v_crimp_actual = m_mix_per_tube / (rho_gas_exit * A_crimp_actual)
            dP_crimp       = rho_gas_exit * v_crimp_actual**2 / (2 * Cd_crimp**2)
        else:
            d_crimp_mm = 0.8
            A_crimp_actual = math.pi * (d_crimp_mm / 2000) ** 2
            v_crimp_actual = m_mix_per_tube / (rho_gas_exit * A_crimp_actual)
            dP_crimp       = rho_gas_exit * v_crimp_actual**2 / (2 * Cd_crimp**2)

        dP_total_Pa    = dP_tube_Pa + dP_crimp
        dP_total_bar   = dP_total_Pa / 1e5
        v_exit_crimped = v_crimp_actual

        self.res['vap_n']                = n_tubes
        self.res['vap_n_raw']            = n_raw
        self.res['vap_od_mm']            = od_tube * 1000
        self.res['vap_id_mm']            = id_tube * 1000
        self.res['vap_crimp_dia_mm']     = round(d_crimp_mm, 2)
        self.res['vap_pitch_actual_mm']  = dome_circumf_m / n_tubes * 1000
        self.res['vap_exit_v_crimped']   = round(v_exit_crimped, 2)
        self.res['vap_bend_radius_mm']   = od_tube * 1000 * 1.5
        self.res['vap_m_air_per_tube']   = m_air_vap / n_tubes * 1000
        self.res['vap_m_fuel_per_tube']  = m_fuel_per_tube * 1000
        self.res['vap_AFR_tube']         = (m_air_vap / n_tubes) / m_fuel_per_tube

        m_air_vap_total  = self.res.get('split_primary_vap', 0.0)
        m_air_vap_per    = m_air_vap_total / n_tubes if n_tubes > 0 else 0.0
        rho_ann          = self.res['rho2']
        v_scoop_target   = self.DESIGN_PARAMS['vap_scoop_target_vel_m_s']
        Cd_scoop         = self.DESIGN_PARAMS['vap_scoop_cd']
        A_scoop_m2       = m_air_vap_per / (rho_ann * v_scoop_target * Cd_scoop)
        d_scoop_mm       = math.sqrt(4 * A_scoop_m2 / math.pi) * 1000
        self.res['vap_scoop_dia_mm']      = round(d_scoop_mm, 2)
        self.res['vap_scoop_vel_actual']  = round(v_scoop_target, 1)

        self.res['fuel_Re_tube']          = Re_tube
        self.res['fuel_dP_tube_Pa']       = dP_tube_Pa
        self.res['fuel_dP_crimp_Pa']      = dP_crimp
        self.res['fuel_dP_total_bar']     = dP_total_bar
        self.res['fuel_dP_target_bar']    = self.DESIGN_PARAMS['fuel_dP_target_bar']
        self.res['vap_tube_flow_regime']  = "Laminar" if Re_tube < 2300 else "Turbulent"
        self.res['vap_mix_velocity']      = round(target_v_mix, 1)

    def pressure_balance_split(self):
        rho2   = self.res['rho2']
        m_air  = self.res['mdot_air']
        m_vap  = self.res.get('split_primary_vap', 0.0)
        mu_air = 1.85e-5

        od_casing = self.inputs['casing_od_inch'] * 0.0254
        wall_m    = self.inputs['wall_thickness_mm'] / 1000.0
        id_casing = od_casing - 2 * wall_m
        od_tunnel = self.inputs['shaft_tunnel_od_inch'] * 0.0254

        r_ol_od  = self.res['outer_liner_od_mm'] / 2 / 1000.0
        r_il_id  = self.res['inner_liner_id_mm'] / 2 / 1000.0
        r_casing_id  = id_casing / 2
        r_tunnel_od  = od_tunnel / 2
        A_outer_geom = math.pi * (r_casing_id**2 - r_ol_od**2)
        A_inner_geom = math.pi * (r_il_id**2     - r_tunnel_od**2)
        D_outer_hyd = id_casing - 2 * r_ol_od
        D_inner_hyd = 2 * (r_il_id - r_tunnel_od)
        L_annulus = self.res['chamber_length_mm'] / 1000.0
        K_turn     = self.DESIGN_PARAMS['K_turn']
        K_entrance = self.DESIGN_PARAMS['K_entrance']

        def darcy_f(Re):
            if Re < 2300:
                return 64.0 / max(Re, 1.0)
            return 0.316 * Re**(-0.25)

        def _path_losses(f_outer):
            m_outer = m_air * f_outer
            m_inner = m_air * (1.0 - f_outer)
            v_outer = m_outer / (rho2 * A_outer_geom) if A_outer_geom > 0 else 0.0
            v_inner = m_inner / (rho2 * A_inner_geom) if A_inner_geom > 0 else 0.0
            Re_outer = rho2 * v_outer * D_outer_hyd / mu_air if D_outer_hyd > 0 else 1.0
            Re_inner = rho2 * v_inner * D_inner_hyd / mu_air if D_inner_hyd > 0 else 1.0
            f_out = darcy_f(Re_outer)
            f_inn = darcy_f(Re_inner)
            dP_fr_out = f_out * (L_annulus / max(D_outer_hyd, 1e-6)) * (0.5 * rho2 * v_outer**2)
            dP_fr_in  = f_inn * (L_annulus / max(D_inner_hyd, 1e-6)) * (0.5 * rho2 * v_inner**2)
            f_extract = m_vap / m_outer if m_outer > 0 else 0.0
            dP_scoop  = 0.5 * rho2 * v_outer**2 * f_extract**2
            dP_turn   = K_turn     * 0.5 * rho2 * v_outer**2
            dP_ent    = K_entrance * 0.5 * rho2 * v_inner**2
            dP_outer = dP_fr_out + dP_scoop + dP_turn
            dP_inner = dP_fr_in  + dP_ent
            return dP_outer, dP_inner

        dP_ref_outer, dP_ref_inner = _path_losses(self.DESIGN_PARAMS['f_outer_feed'])
        dP_ref_mean  = (dP_ref_outer + dP_ref_inner) / 2.0
        tol          = max(5.0, dP_ref_mean * 0.01)

        # Physics floor: outer path must carry at least the vaporizer scoop demand.
        # Previously hard-coded to max(scoop_frac*1.10, 0.45), which overconstrained the
        # search range and hid the true zero-crossing at f≈0.38 for the KJ66 geometry.
        # Correct floor: scoop fraction × 1.10 safety margin, minimum 0.30 (structural).
        # The 0.45 floor was conservative and physically unjustified.
        f_lo_physics = (m_vap / m_air) * 1.10
        f_lo = max(f_lo_physics, 0.30)
        f_hi = 0.95
        converged  = False
        f_solved   = self.DESIGN_PARAMS['f_outer_feed']

        iter_count = 0
        try:
            def residual(f):
                o, i = _path_losses(f)
                return o - i
            r_lo = residual(f_lo)
            r_hi = residual(f_hi)
            if r_lo * r_hi < 0:
                # True zero-crossing found — use bisection
                for _iter in range(50):
                    iter_count = _iter + 1
                    f_mid = (f_lo + f_hi) / 2.0
                    r_mid = residual(f_mid)
                    if abs(r_mid) < tol:
                        f_solved  = f_mid
                        converged = True
                        break
                    if r_lo * r_mid < 0:
                        f_hi, r_hi = f_mid, r_mid
                    else:
                        f_lo, r_lo = f_mid, r_mid
                else:
                    f_solved  = (f_lo + f_hi) / 2.0
                    converged = abs(residual(f_solved)) < tol * 5.0
            else:
                # No zero-crossing in valid range — use golden-section search to find
                # the split that minimises |ΔP_outer − ΔP_inner| (closest to balance).
                # This is physically the best-achievable split given the geometry.
                phi = (math.sqrt(5) - 1) / 2  # golden ratio
                a, b = f_lo, f_hi
                for _iter in range(60):
                    iter_count = _iter + 1
                    c = b - phi * (b - a)
                    d = a + phi * (b - a)
                    if abs(residual(c)) < abs(residual(d)):
                        b = d
                    else:
                        a = c
                    if (b - a) < 1e-6:
                        break
                f_solved  = (a + b) / 2.0
                converged = False   # flag as fallback — not a true pressure balance
        except Exception as e:
            self.res['split_solver_error'] = str(e)

        dP_outer_final, dP_inner_final = _path_losses(f_solved)
        final_residual = abs(dP_outer_final - dP_inner_final)
        self.res['f_outer_solved']         = f_solved
        self.res['split_converged']        = converged
        self.res['dP_outer_path_Pa']       = dP_outer_final
        self.res['dP_inner_path_Pa']       = dP_inner_final
        self.res['split_solver_tol_Pa']    = round(tol, 1)
        self.res['split_solver_iters']     = iter_count
        self.res['split_solver_residual_Pa'] = round(final_residual, 1)

    def combustion_loading(self):
        P3_kPa  = self.res['P2_Pa'] / 1000.0
        T3      = self.res['T2_K']
        m_air   = self.res['mdot_air']
        A_comb    = self.res['combustion_annulus_A']
        L_primary = self.res['L_primary_mm'] / 1000.0
        V_primary = A_comb * L_primary
        V_total   = A_comb * self.res['chamber_length_mm'] / 1000.0
        CLP   = m_air / (P3_kPa * V_primary)
        CLP_T = CLP * math.sqrt(T3 / 300.0)

        # Combustor residence time τ = V_total / V̇
        # V̇ = ṁ_total / ρ_mean  using the mean combustion zone density.
        # Mean temperature = (T2 + TIT) / 2, which gives a realistic bulk average
        # across the primary and secondary zones. This is the correct reference
        # density for estimating how long the mixture spends in the chamber.
        #
        # Note: using cold inlet density (ρ_inlet) over-estimates τ by ~2× because
        # cold air occupies less volume than hot combustion gases — the comment in
        # the previous version labelling this "conservative" was incorrect.
        #
        # Industry target: τ > 2 ms at mean conditions for stable combustion at
        # micro-GT pressures. Flag below 1.5 ms: high blowout risk at part throttle.
        T_mean_comb = (self.res['T2_K'] + self.inputs['target_tit_k']) / 2.0
        rho_mean    = self.res['P2_Pa'] / (self.R * T_mean_comb)
        m_total_comb = m_air + self.res['mdot_fuel']
        V_dot       = m_total_comb / rho_mean
        tau_comb_ms = (V_total / V_dot) * 1000.0
        self.res['CLP']              = CLP
        self.res['CLP_T']            = CLP_T
        self.res['CLP_P3_kPa']       = P3_kPa
        self.res['CLP_V_primary_m3'] = V_primary
        self.res['CLP_V_total_m3']   = V_total
        self.res['tau_comb_ms']      = round(tau_comb_ms, 2)
        self.res['CLP_stable']       = CLP < 8.0
        self.res['CLP_acceptable']   = CLP < 15.0

        # Combustion efficiency estimate from CLP (Lefebvre, GTC 3rd Ed. Table 5.1)
        # At CLP < 10: η > 0.999; CLP = 15: η ≈ 0.97; CLP = 20: η ≈ 0.94; CLP = 25: η ≈ 0.90
        # Linear fit: η = max(0.85, 1.0 - 0.006*(CLP-10)) for CLP > 10
        # This aligns with the fixed 0.96 used in mass_flow_and_fuel (CLP≈15 → η=0.97)
        eta_comb_estimated = max(0.85, min(0.999, 1.0 - 0.006 * max(0.0, CLP - 10.0)))
        self.res['eta_comb_estimated'] = round(eta_comb_estimated, 3)
    def exit_conditions(self):
        """
        1-D turbine inlet (combustor exit) conditions.
        Uses γ=1.33 for hot combustion products instead of the cold-air 1.4
        used elsewhere in the code — matters for Mach and isentropic relations.
        P4 is derived from the liner ΔP design target (not independently verified).
        """
        GAMMA_HOT = 1.33          # combustion products at ~900–1200 K
        R         = self.R        # 287.05 J/kg·K — close enough for lean products

        P2        = self.res['P2_Pa']
        T4        = self.inputs['target_tit_k']
        dP_frac   = self.DESIGN_PARAMS['target_pressure_drop']
        P4        = P2 * (1.0 - dP_frac)

        mdot_exit = self.res['mdot_air'] + self.res['mdot_fuel']
        A_exit    = self.res['combustion_annulus_A']
        cp4       = self._get_cp(T4)

        rho4  = P4 / (R * T4)
        V4    = mdot_exit / (rho4 * A_exit)
        a4    = math.sqrt(GAMMA_HOT * R * T4)
        Ma4   = V4 / a4

        # Isentropic total conditions (stagnation)
        factor        = 1.0 + (GAMMA_HOT - 1.0) / 2.0 * Ma4 ** 2
        T4_total      = T4 * factor
        P4_total      = P4 * factor ** (GAMMA_HOT / (GAMMA_HOT - 1.0))

        # Specific enthalpy available to the turbine (total enthalpy drop to ambient)
        T_amb    = 288.15
        cp_mean  = self._get_cp((T4 + T_amb) / 2.0)
        h_avail  = cp_mean * (T4_total - T_amb)   # J/kg — upper bound on turbine work

        self.res['exit_P4_Pa']       = P4
        self.res['exit_P4_kPa']      = P4 / 1000.0
        self.res['exit_T4_K']        = T4
        self.res['exit_rho4']        = rho4
        self.res['exit_V4_m_s']      = V4
        self.res['exit_Ma4']         = Ma4
        self.res['exit_a4_m_s']      = a4
        self.res['exit_T4_total_K']  = T4_total
        self.res['exit_P4_total_Pa'] = P4_total
        self.res['exit_mdot_kg_s']   = mdot_exit
        self.res['exit_cp4']         = cp4
        self.res['exit_h_avail_kJ_kg'] = h_avail / 1000.0
        self.res['exit_gamma_hot']   = GAMMA_HOT
        self.res['exit_Ma4_high']    = Ma4 > 0.25   # flag: compressibility matters above ~0.25

    def hole_sizing(self):
        dP   = self.res['P2_Pa'] * self.DESIGN_PARAMS['target_pressure_drop']
        Cd   = self.DESIGN_PARAMS['discharge_coeff_hole']
        rho2 = self.res['rho2']
        flow_factor = Cd * math.sqrt(2 * rho2 * dP)

        A_pri = self.res['split_primary_liner'] / flow_factor
        A_sec = self.res['split_secondary']     / flow_factor
        A_dil = self.res['split_dilution']      / flow_factor

        if self.DESIGN_PARAMS['outer_liner_hole_fraction'] is not None:
            f_outer = self.DESIGN_PARAMS['outer_liner_hole_fraction']
        else:
            C_outer = self.res['outer_liner_circumf_m']
            C_inner = self.res['inner_liner_circumf_m']
            f_outer = C_outer / (C_outer + C_inner)
        f_inner = 1.0 - f_outer

        n_vap = self.res['vap_n']

        def hole_dia(area_total, fraction, n_holes):
            A_each = (area_total * fraction) / n_holes
            return math.sqrt(4 * A_each / math.pi)

        n_pri_out = n_vap * 2
        n_pri_in  = n_vap * 2

        sec_mult  = max(1, int(self.DESIGN_PARAMS['sec_holes_per_vap']))
        n_sec_out = n_vap * sec_mult
        n_sec_in  = n_vap * sec_mult

        # Dilution: configurable multiplier via dil_holes_per_vap (default 4×/vap).
        # At 4×/vap (24 holes/liner at n_vap=6): d_dil_out ≈ 5.6 mm, penetration ~1.09×H.
        # Combined coverage ≈ 1.97×H — full annulus swept by opposing jets.
        # KJ66 historical match: set dil_holes_per_vap=2 in the KJ66 preset.
        dil_mult  = max(1, int(self.DESIGN_PARAMS.get('dil_holes_per_vap', 4)))
        n_dil_out = n_vap * dil_mult
        n_dil_in  = n_vap * dil_mult

        d_pri_out = hole_dia(A_pri, f_outer, n_pri_out)
        d_pri_in  = hole_dia(A_pri, f_inner, n_pri_in)
        d_sec_out = hole_dia(A_sec, f_outer, n_sec_out)
        d_sec_in  = hole_dia(A_sec, f_inner, n_sec_in)
        d_dil_out = hole_dia(A_dil, f_outer, n_dil_out)
        d_dil_in  = hole_dia(A_dil, f_inner, n_dil_in)

        T2      = self.res['T2_K']
        a_sound = math.sqrt(self.GAMMA * self.R * T2)

        def hole_mach(m_zone, n_holes, d_hole):
            A_hole = math.pi * (d_hole / 2) ** 2
            if A_hole <= 0 or n_holes <= 0:
                return 0.0
            v_hole = m_zone / (rho2 * A_hole * n_holes)
            return v_hole / a_sound

        def compressible_correction(Ma):
            return math.sqrt(1.0 + 0.2 * Ma**2)

        Ma_pri_0 = hole_mach(self.res['split_primary_liner'], n_pri_out, d_pri_out)
        Ma_sec_0 = hole_mach(self.res['split_secondary'],     n_sec_out, d_sec_out)
        Ma_dil_0 = hole_mach(self.res['split_dilution'],      n_dil_out, d_dil_out)

        def corrected_dia(d_incomp, Ma):
            factor = compressible_correction(Ma)
            A_corr = math.pi * (d_incomp / 2)**2 * factor
            return math.sqrt(4 * A_corr / math.pi)

        d_pri_out = corrected_dia(d_pri_out, Ma_pri_0)
        d_pri_in  = corrected_dia(d_pri_in,  Ma_pri_0)
        d_sec_out = corrected_dia(d_sec_out, Ma_sec_0)
        d_sec_in  = corrected_dia(d_sec_in,  Ma_sec_0)
        d_dil_out = corrected_dia(d_dil_out, Ma_dil_0)
        d_dil_in  = corrected_dia(d_dil_in,  Ma_dil_0)

        Ma_pri = hole_mach(self.res['split_primary_liner'], n_pri_out, d_pri_out)
        Ma_sec = hole_mach(self.res['split_secondary'],     n_sec_out, d_sec_out)
        Ma_dil = hole_mach(self.res['split_dilution'],      n_dil_out, d_dil_out)

        self.res['Ma_hole_pri'] = Ma_pri
        self.res['Ma_hole_sec'] = Ma_sec
        self.res['Ma_hole_dil'] = Ma_dil
        self.res['Ma_hole_pri_uncorr'] = Ma_pri_0
        self.res['Ma_hole_sec_uncorr'] = Ma_sec_0
        self.res['Ma_hole_dil_uncorr'] = Ma_dil_0
        self.res['hole_compressible_warning'] = any(m > 0.3 for m in [Ma_pri_0, Ma_sec_0, Ma_dil_0])

        m_film       = self.res['m_film_cooling']
        A_film_total = m_film / flow_factor
        n_rows_pri = self.DESIGN_PARAMS['film_rows_primary']
        n_rows_sec = self.DESIGN_PARAMS['film_rows_secondary']
        n_rows_dil = self.DESIGN_PARAMS['film_rows_dilution']
        n_rows     = n_rows_pri + n_rows_sec + n_rows_dil
        A_per_row   = A_film_total / max(n_rows, 1)
        d_film_m    = self.DESIGN_PARAMS['film_hole_dia_mm'] / 1000.0
        a_film_hole = math.pi * (d_film_m / 2) ** 2
        n_film_out_per_row = math.ceil(A_per_row * f_outer / a_film_hole)
        n_film_in_per_row  = math.ceil(A_per_row * f_inner / a_film_hole)
        pitch_m = self.DESIGN_PARAMS['film_row_pitch_mm'] / 1000.0
        row_positions = [pitch_m * (i + 0.5) for i in range(n_rows)]

        self.res['film_n_rows']              = n_rows
        self.res['film_n_rows_pri']          = n_rows_pri
        self.res['film_n_rows_sec']          = n_rows_sec
        self.res['film_n_rows_dil']          = n_rows_dil
        self.res['film_hole_dia_mm_actual']  = self.DESIGN_PARAMS['film_hole_dia_mm']
        self.res['film_holes_per_row_outer'] = n_film_out_per_row
        self.res['film_holes_per_row_inner'] = n_film_in_per_row
        self.res['film_total_area_mm2']      = round(A_film_total * 1e6, 1)
        self.res['film_total_holes']         = n_rows * (n_film_out_per_row + n_film_in_per_row)

        T_gas_dil  = 1200.0
        rho_gas    = self.res['P2_Pa'] / (self.R * T_gas_dil)
        m_at_dil   = self.res['split_primary'] + self.res['split_secondary'] + self.res['mdot_fuel']
        A_comb     = self.res['combustion_annulus_A']
        cp_dil     = self._get_cp(T_gas_dil)
        V_gas      = m_at_dil / (rho_gas * A_comb)

        A_hole_out = math.pi * (d_dil_out / 2) ** 2
        V_jet_out  = (self.res['split_dilution'] * f_outer / n_dil_out) / (rho2 * A_hole_out)
        A_hole_in  = math.pi * (d_dil_in / 2) ** 2
        V_jet_in   = (self.res['split_dilution'] * f_inner / n_dil_in) / (rho2 * A_hole_in)

        J_out = (rho2 * V_jet_out ** 2) / (rho_gas * V_gas ** 2) if V_gas > 0 else 0
        J_in  = (rho2 * V_jet_in  ** 2) / (rho_gas * V_gas ** 2) if V_gas > 0 else 0

        H_comb = self.res['combustion_gap_mm'] / 1000
        pen_out_norm = 1.15 * math.sqrt(J_out) * (d_dil_out / H_comb) if H_comb > 0 else 0
        pen_in_norm  = 1.15 * math.sqrt(J_in)  * (d_dil_in  / H_comb) if H_comb > 0 else 0
        pen_out_mm   = pen_out_norm * H_comb * 1000
        pen_in_mm    = pen_in_norm  * H_comb * 1000

        self.res['hole_split_f_outer']    = f_outer
        self.res['hole_split_f_inner']    = f_inner
        self.res['pri_out_qty']  = n_pri_out
        self.res['pri_out_mm']   = d_pri_out * 1000
        self.res['sec_out_qty']  = n_sec_out
        self.res['sec_out_mm']   = d_sec_out * 1000
        self.res['dil_out_qty']  = n_dil_out
        self.res['dil_out_mm']   = d_dil_out * 1000
        self.res['pri_in_qty']   = n_pri_in
        self.res['pri_in_mm']    = d_pri_in  * 1000
        self.res['sec_in_qty']   = n_sec_in
        self.res['sec_in_mm']    = d_sec_in  * 1000
        self.res['dil_in_qty']   = n_dil_in
        self.res['dil_in_mm']    = d_dil_in  * 1000
        self.res['dil_J_outer']         = J_out
        self.res['dil_J_inner']         = J_in
        self.res['dil_pen_outer_mm']    = pen_out_mm
        self.res['dil_pen_inner_mm']    = pen_in_mm
        self.res['H_comb_mm']           = H_comb * 1000
        self.res['dil_pen_norm_outer']  = pen_out_norm
        self.res['dil_pen_norm_inner']  = pen_in_norm

    def liner_structural(self):
        wall_m = self.inputs['wall_thickness_mm'] / 1000.0
        dP     = self.res['P2_Pa'] * self.DESIGN_PARAMS['target_pressure_drop']
        r_outer_mean = (self.res['outer_liner_od_mm'] + self.res['outer_liner_id_mm']) / 2 / 1000 / 2
        sigma_outer  = dP * r_outer_mean / wall_m / 1e6
        r_inner_mean = (self.res['inner_liner_od_mm'] + self.res['inner_liner_id_mm']) / 2 / 1000 / 2
        sigma_inner  = dP * r_inner_mean / wall_m / 1e6
        sigma_allow_ss304   =  65.0
        sigma_allow_in625   = 175.0
        circ_outer_mm = self.res['outer_liner_circumf_m'] * 1000
        pitch_out_mm  = circ_outer_mm / self.res['dil_out_qty']
        lig_out_mm    = pitch_out_mm - self.res['dil_out_mm']
        circ_inner_mm = self.res['inner_liner_circumf_m'] * 1000
        pitch_in_mm   = circ_inner_mm / self.res['dil_in_qty']
        lig_in_mm     = pitch_in_mm - self.res['dil_in_mm']
        lig_min_req   = max(self.inputs['wall_thickness_mm'] * 2.0, 1.5)
        self.res['sigma_hoop_outer_MPa']  = sigma_outer
        self.res['sigma_hoop_inner_MPa']  = sigma_inner
        self.res['sigma_allow_ss304_MPa'] = sigma_allow_ss304
        self.res['sigma_allow_in625_MPa'] = sigma_allow_in625
        self.res['ligament_outer_mm']     = lig_out_mm
        self.res['ligament_inner_mm']     = lig_in_mm
        self.res['ligament_min_req_mm']   = lig_min_req
        self.res['liner_ss304_ok']        = (sigma_outer < sigma_allow_ss304 and sigma_inner < sigma_allow_ss304)
        self.res['liner_in625_ok']        = (sigma_outer < sigma_allow_in625 and sigma_inner < sigma_allow_in625)
        self.res['ligament_outer_ok']     = lig_out_mm >= lig_min_req
        self.res['ligament_inner_ok']     = lig_in_mm  >= lig_min_req

    def temperature_traverse_quality(self):
        pen_norm_out = self.res['dil_pen_norm_outer']
        pen_norm_in  = self.res['dil_pen_norm_inner']
        eff_out = 1.0 - abs(pen_norm_out - 0.5) / 0.5
        eff_in  = 1.0 - abs(pen_norm_in  - 0.5) / 0.5
        eff_out = max(0.0, min(1.0, eff_out))
        eff_in  = max(0.0, min(1.0, eff_in))
        mix_eff = (eff_out + eff_in) / 2.0
        MI = 0.40 - mix_eff * 0.35
        combined_coverage = pen_norm_out + pen_norm_in
        self.res['mixing_index']       = MI
        self.res['mix_eff_outer']      = eff_out
        self.res['mix_eff_inner']      = eff_in
        self.res['dil_combined_cov']   = combined_coverage

    def get_cad_geometry(self):
        return {
            "metadata": {
                "units": "mm",
                "description": "Reverse-flow annular micro-jet combustor — KJ66 style",
                "generated_by": "V21_CombustionChamberDesign.py",
                "liner_material": self.res.get('liner_material', '304SS'),
                "liner_material_name": self.res.get('liner_material_name', '304 Stainless Steel'),
                "liner_alpha_per_K": self.res.get('liner_alpha_per_K', 17.2e-6),
                "heat_loss_factor": self.res.get('heat_loss_factor', 1.00),
                "dimension_note": (
                    "Liner OD values are COLD BUILD dimensions (subtract thermal offset). "
                    "Parts sized at room temperature will expand to nominal hot dimensions "
                    "at operating temperature (~900 C). Do NOT use od_hot in CAD for fabrication."
                )
            },
            "overall": {"chamber_length": self.res['chamber_length_mm'], "mean_comb_dia": self.res['D_mean_comb_mm']},
            "casing": {"od": self.res['casing_od_mm'], "id": self.res['casing_id_mm'], "wall_thickness": self.inputs['wall_thickness_mm']},
            "shaft_tunnel": {"od": self.res['shaft_tunnel_od_mm'], "id": self.res['shaft_tunnel_id_mm']},
            "outer_liner": {"od": self.res['outer_liner_od_cold_mm'], "od_hot": self.res['outer_liner_od_hot_mm'], "thermal_offset": self.res['outer_liner_thermal_offset_mm'], "id": self.res['outer_liner_id_mm'], "length": self.res['chamber_length_mm'], "annulus_gap": self.res['outer_annulus_gap_mm']},
            "inner_liner": {"od": self.res['inner_liner_od_cold_mm'], "od_hot": self.res['inner_liner_od_hot_mm'], "thermal_offset": self.res['inner_liner_thermal_offset_mm'], "id": self.res['inner_liner_id_mm'], "length": self.res['chamber_length_mm'], "annulus_gap": self.res['inner_annulus_gap_mm']},
            "combustion_annulus": {"radial_height": self.res['combustion_gap_mm'], "mean_diameter": self.res['D_mean_comb_mm'], "cross_section_area_mm2": self.res['combustion_annulus_A'] * 1e6},
            "vaporizers": {"count": self.res['vap_n'], "raw_count": self.res['vap_n_raw'], "pitch_actual": self.res['vap_pitch_actual_mm'], "tube_od": self.res['vap_od_mm'], "tube_bore_id": self.res['vap_id_mm'], "crimp_orifice_id": self.res['vap_crimp_dia_mm'], "scoop_inlet_dia": self.res['vap_scoop_dia_mm'], "bend_radius": self.res['vap_bend_radius_mm'], "per_tube": {"fuel_g_s": self.res['vap_m_fuel_per_tube'], "air_g_s": self.res['vap_m_air_per_tube']}},
            "main_holes": {"outer_fraction": self.res['hole_split_f_outer'], "inner_fraction": self.res['hole_split_f_inner'], "primary": {"outer_qty": self.res['pri_out_qty'], "outer_dia": self.res['pri_out_mm'], "inner_qty": self.res['pri_in_qty'], "inner_dia": self.res['pri_in_mm']}, "secondary": {"outer_qty": self.res['sec_out_qty'], "outer_dia": self.res['sec_out_mm'], "inner_qty": self.res['sec_in_qty'], "inner_dia": self.res['sec_in_mm']}, "dilution": {"outer_qty": self.res['dil_out_qty'], "outer_dia": self.res['dil_out_mm'], "inner_qty": self.res['dil_in_qty'], "inner_dia": self.res['dil_in_mm']}},
            "film_cooling": {"total_rows": self.res['film_n_rows'], "rows_primary": self.res['film_n_rows_pri'], "rows_secondary": self.res['film_n_rows_sec'], "rows_dilution": self.res['film_n_rows_dil'], "row_pitch": self.DESIGN_PARAMS['film_row_pitch_mm'], "hole_dia": self.res['film_hole_dia_mm_actual'], "outer_holes_per_row": self.res['film_holes_per_row_outer'], "inner_holes_per_row": self.res['film_holes_per_row_inner'], "total_holes": self.res['film_total_holes'], "total_area_mm2": self.res['film_total_area_mm2']},
            "annulus_velocities": {"outer_entry": self.res['v_outer_annulus'], "outer_post_vap": self.res['v_outer_post_vap'], "inner": self.res['v_inner_annulus']},
            "zone_lengths": {"primary": self.res['L_primary_mm'], "secondary": self.res['L_secondary_mm'], "dilution": self.res['L_dilution_mm']}
        }

    def stability_checks(self):
        """
        Three fluid-dynamic stability checks that a 1-D sizing code cannot guarantee
        but CAN flag. These are the failure modes that kill most reverse-flow microjets
        even when the static geometry looks correct.

        CHECK 1 — PRIMARY JET RECIRCULATION STRENGTH (J_primary)
        ──────────────────────────────────────────────────────────
        The toroidal recirculation vortex in the primary zone is maintained by the
        radial momentum of the primary jets opposing the axial bulk flow.

        J = (rho_jet * V_jet²) / (rho_hot * V_axial²)
          rho_jet  = compressor exit density (cold feed air entering hole)
          V_jet    = primary liner hole jet velocity
          rho_hot  = primary zone gas density at estimated T_primary
          V_axial  = bulk axial velocity through combustion annulus cross-section

        Thresholds (Lefebvre §4.3, adjusted for dual-wall reverse-flow geometry):
          J < 5   : vortex collapse risk — jets swept downstream, flame lifts off
          J 5–80  : stable recirculation (note: dual-wall and vaporizer assist lower bound)
          J > 80  : jets over-penetrate, vortex splits into two counter-rotating cells

        CHECK 2 — VAPORIZER FLASHBACK VELOCITY
        ───────────────────────────────────────
        Flashback occurs when the fuel/air mixture velocity inside the tube falls below
        the laminar flame speed of the mixture. For Jet-A/kerosene air at vaporizer
        conditions (phi~0.8–1.2, T~450–700K):
          S_L ≈ 0.35–0.6 m/s (Metgalchi & Keck 1982, kerosene surrogate)
        Safe design: V_mixture > 10× S_L ≈ 4 m/s (Lefebvre §10.2)
        Critical:    V_mixture < 2× S_L ≈ 0.8 m/s (immediate flashback risk)

        Also check tube L/D: if L/D < 8 and D > 6mm, quench distance is insufficient.

        CHECK 3 — DOME HEAT FLUX ESTIMATE
        ──────────────────────────────────
        The dome sits in the primary recirculation bubble and sees ~T_primary with
        minimal cooling air. Heat flux via forced convection from impinging jets:
          Nu = 0.037 * Re^0.8 * Pr^0.33  (turbulent impingement, Dittus-Boelter variant)
          h  = Nu * k_gas / d_jet
          q  = h * (T_gas - T_wall)

        Limits (uninsulated 304SS or similar):
          q < 300 kW/m²  : acceptable, passive radiation + film cooling manages
          300–600 kW/m²  : high risk — add dome cooling air distribution ring
          > 600 kW/m²    : critical — dome will fail without dedicated transpiration cooling
        """
        # ── 1. J_primary ──────────────────────────────────────────────────────
        rho2    = self.res['rho2']
        T_pri   = self.res.get('T_primary_zone_est', 1800.0)
        P2      = self.res['P2_Pa']
        rho_pri = P2 / (self.R * T_pri)

        # Primary jet velocity: liner-hole fraction of primary air through outer holes
        m_air_pri_liner = self.res['split_primary_liner']
        n_pri_out       = self.res['pri_out_qty']
        d_pri_out_m     = self.res['pri_out_mm'] / 1000.0
        A_pri_out       = math.pi * (d_pri_out_m / 2) ** 2
        V_jet_pri       = (m_air_pri_liner * self.res['hole_split_f_outer']
                           / (rho2 * A_pri_out * n_pri_out)) if A_pri_out > 0 else 0.0

        # Axial bulk velocity in combustion annulus (cold density — conservative)
        A_comb  = self.res['combustion_annulus_A']
        m_comb  = self.res['m_combustion_air']
        V_axial = m_comb / (rho2 * A_comb) if A_comb > 0 else 1.0

        J_primary = (rho2 * V_jet_pri ** 2) / (rho_pri * V_axial ** 2) if V_axial > 0 else 0.0

        # ── 2. Vaporizer flashback ─────────────────────────────────────────────
        n_tubes       = self.res['vap_n']
        m_fuel_total  = self.res['mdot_fuel']
        m_air_vap     = self.res['split_primary_vap']
        m_mix_tube    = (m_fuel_total + m_air_vap) / n_tubes  # kg/s per tube

        id_tube_m     = self.res['vap_id_mm'] / 1000.0
        A_tube        = math.pi * (id_tube_m / 2) ** 2

        # Mixture density inside vaporizer at ~vaporization temperature
        T_vap_zone    = 600.0      # K — conservative estimate (fuel + air mixing zone)
        rho_mix       = P2 / (self.R * T_vap_zone)
        V_mix         = m_mix_tube / (rho_mix * A_tube) if A_tube > 0 else 0.0

        # L/D quench check
        L_tube_m = (self.res.get('chamber_length_mm', 120) / 1000) * 1.5 + 2 * (self.res['vap_od_mm'] / 1000)
        LD_tube  = L_tube_m / id_tube_m if id_tube_m > 0 else 0.0

        # Vaporizer residence time: τ_vap = L_tube / V_mix
        # Reference threshold for air-assisted tube vaporizers:
        # With co-flowing hot air and a tube wall pre-heated by recirculation gases,
        # Jet-A vaporizes in <0.5 ms at these conditions (Colket & Spadaccini 2001).
        # The 2–5 ms reference in earlier versions applied to UNASSISTED liquid spray
        # combustors — not applicable here. Conservative threshold: 0.3 ms.
        # Below 0.3 ms: tube is too short for adequate fuel/air mixing before the
        # 180° bend; partially-mixed charge risks stratified combustion and coking.
        tau_vap_ms = (L_tube_m / V_mix * 1000.0) if V_mix > 0 else 0.0

        # ── 3. Dome heat flux — convective + radiative ───────────────────────
        T_wall_dome   = self.inputs['target_tit_k'] * 0.80   # dome sheet ~80% TIT (conservative)
        T_gas         = max(T_pri, 1200.0)

        mu_gas = 3.5e-5 * (T_gas / 1000.0) ** 0.7
        k_gas  = 0.055 * (T_gas / 1000.0) ** 0.8
        Pr     = 0.72
        d_jet_m = d_pri_out_m
        Re_imp  = rho2 * V_jet_pri * d_jet_m / mu_gas if mu_gas > 0 else 0.0
        # Physical basis: the impinging jet is cold feed air exiting the primary liner holes
        # from the outer feed annulus at density rho2. Using rho_pri (hot primary zone gas,
        # ~5× less dense) understates Re by 5×, Nu by 3.7×, and q_conv by 3.7×.
        # The jet heats AFTER impingement — at the hole exit it is still at inlet conditions.
        Nu_imp  = 0.037 * Re_imp ** 0.8 * Pr ** 0.33 if Re_imp > 0 else 0.0
        h_conv  = Nu_imp * k_gas / d_jet_m if d_jet_m > 0 else 0.0
        q_conv  = h_conv * (T_gas - T_wall_dome) / 1000.0   # kW/m²

        # Radiative component — significant in rich primary zones.
        # Emissivity scales with primary zone richness: more soot at higher phi.
        # phi=1.0 (stoich) → ε=0.30; phi=1.6 (KJ66 rich) → ε=0.36; clamped 0.30–0.50.
        # Ref: Lefebvre §9.3 — soot formation increases optical depth above phi~1.2.
        phi_pri_act   = self.res.get('phi_primary_actual', 1.6)
        epsilon_flame = max(0.30, min(0.50, 0.30 + 0.10 * (phi_pri_act - 1.0)))

        # View factor: the flame does not completely envelop the dome (F=1.0 assumption).
        # In a reverse-flow annular combustor: vaporizer tubes, dome lip, and cooling gap
        # reduce the effective radiation view. F_view ≈ 0.4–0.7; 0.60 is a conservative
        # midpoint for an unshielded dome. Without this, q_rad is overpredicted by ~67%.
        # Ref: Lefebvre §9.3 — typical annular combustor dome view factor range.
        F_view   = 0.60
        sigma_sb = 5.6704e-8   # W/m²K⁴
        q_rad    = F_view * epsilon_flame * sigma_sb * (T_gas**4 - T_wall_dome**4) / 1000.0  # kW/m²
        q_dome   = q_conv + q_rad

        # ── 4. Vaporizer tube thermal check ──────────────────────────────────
        # 1-D thermal resistance model: flame→tube_wall→mixture
        # Outer surface: recirculation zone impingement, h_outer ≈ 300 W/m²K
        # Inner surface: turbulent mixture flow, h_inner from Dittus-Boelter
        # This determines whether the tube material choice is safe and flags coking risk.
        #
        # Jet-A coking threshold: ~650 K tube wall temperature (inner surface)
        # Above this, aromatic hydrocarbons polymerise → carbon deposits → tube blockage
        # KJ66 builders report vaporizer replacement every 5–20 flight hours — this IS
        # the dominant life-limiting mechanism. The model quantifies the risk.
        T_mix_inside = 600.0   # K — conservative mixture temperature inside tube
        h_outer_vap  = 300.0   # W/m²K — recirculation zone impingement (conservative)

        # Inner h from mixture flow (Dittus-Boelter, turbulent heating)
        rho_mix_vap = P2 / (self.R * T_mix_inside)
        mu_mix      = 3.0e-5   # Pa·s — air at 600K
        k_mix       = 0.046    # W/mK — air at 600K
        Re_vap      = rho_mix_vap * V_mix * id_tube_m / mu_mix if (mu_mix > 0 and id_tube_m > 0) else 0.0
        Pr_vap      = 0.71
        Nu_vap      = 0.023 * max(Re_vap, 1.0)**0.8 * Pr_vap**0.4 if Re_vap > 10 else 3.66
        h_inner_vap = Nu_vap * k_mix / id_tube_m if id_tube_m > 0 else 0.0

        t_wall_vap  = 0.0006   # m — 0.6 mm wall (typical 1/5" SS tube)
        k_wall_vap  = 16.0     # W/mK — 316SS / Inconel 625 (similar)
        R_vap       = 1.0/h_outer_vap + t_wall_vap/k_wall_vap + (1.0/h_inner_vap if h_inner_vap > 0 else 0.1)
        q_vap_net   = (T_pri - T_mix_inside) / R_vap if R_vap > 0 else 0.0
        T_vap_wall_outer = T_pri - q_vap_net / h_outer_vap
        T_vap_wall_inner = T_mix_inside + q_vap_net / h_inner_vap if h_inner_vap > 0 else T_vap_wall_outer

        # ── Store all results ──────────────────────────────────────────────────
        self.res['stab_J_primary']           = J_primary
        self.res['stab_V_jet_pri']           = V_jet_pri
        self.res['stab_V_axial']             = V_axial
        self.res['stab_rho_pri']             = rho_pri
        self.res['stab_J_ok']               = 5.0 <= J_primary <= 80.0
        self.res['stab_J_collapse']         = J_primary < 5.0
        self.res['stab_J_overpenetrate']    = J_primary > 80.0
        self.res['stab_V_axial_high']       = V_axial > 30.0

        self.res['stab_V_mix']              = V_mix
        self.res['stab_tau_vap_ms']         = round(tau_vap_ms, 2)
        self.res['stab_LD_tube']            = LD_tube
        self.res['stab_flashback_critical'] = V_mix < 0.8
        self.res['stab_flashback_risk']     = V_mix < 4.0
        self.res['stab_LD_ok']             = LD_tube >= 8.0 or id_tube_m <= 0.006

        self.res['stab_q_conv_kW_m2']      = round(q_conv, 1)
        self.res['stab_q_rad_kW_m2']       = round(q_rad,  1)
        self.res['stab_q_dome_kW_m2']      = round(q_dome, 1)
        self.res['stab_epsilon_flame']     = round(epsilon_flame, 3)
        self.res['stab_F_view']            = F_view
        self.res['stab_T_wall_dome_K']     = T_wall_dome
        self.res['stab_Re_imp']            = Re_imp
        self.res['stab_h_conv']            = h_conv
        self.res['stab_dome_ok']           = q_dome < 300.0
        self.res['stab_dome_high']         = 300.0 <= q_dome < 600.0
        self.res['stab_dome_critical']     = q_dome >= 600.0

        self.res['stab_T_vap_wall_outer_K'] = round(T_vap_wall_outer, 0)
        self.res['stab_T_vap_wall_inner_K'] = round(T_vap_wall_inner, 0)
        self.res['stab_h_inner_vap']        = round(h_inner_vap, 0)
        self.res['stab_vap_coking_risk']    = T_vap_wall_inner > 650.0
        self.res['stab_vap_structural_risk']= T_vap_wall_outer > 1100.0

    def run(self):
        self.thermodynamics()
        self.mass_flow_and_fuel()
        self.zonal_analysis()
        self.mechanical_geometry()
        self.vaporizer_tubes()
        self.pressure_balance_split()
        self._update_annulus_velocities_post_solver()
        self.hole_sizing()
        self.liner_structural()
        self.temperature_traverse_quality()
        self.combustion_loading()
        self.exit_conditions()
        self.stability_checks()
        self.res['cad_geometry'] = self.get_cad_geometry()
        return self.res

    def _update_annulus_velocities_post_solver(self):
        f_solved  = self.res.get('f_outer_solved', self.DESIGN_PARAMS['f_outer_feed'])
        rho2      = self.res['rho2']
        m_air     = self.res['mdot_air']
        m_vap     = self.res.get('split_primary_vap', 0.0)
        wall_m    = self.inputs['wall_thickness_mm'] / 1000.0
        od_casing = self.inputs['casing_od_inch'] * 0.0254
        id_casing = od_casing - 2 * wall_m
        od_tunnel = self.inputs['shaft_tunnel_od_inch'] * 0.0254
        r_casing_id  = id_casing / 2
        r_ol_od      = self.res['outer_liner_od_mm'] / 2 / 1000.0
        r_il_id      = self.res['inner_liner_id_mm'] / 2 / 1000.0
        r_tunnel_od  = od_tunnel / 2
        A_outer = math.pi * (r_casing_id**2 - r_ol_od**2)
        A_inner = math.pi * (r_il_id**2     - r_tunnel_od**2)
        m_outer = m_air * f_solved
        m_inner = m_air * (1.0 - f_solved)
        v_outer       = m_outer / (rho2 * A_outer) if A_outer > 0 else 0.0
        v_inner       = m_inner / (rho2 * A_inner) if A_inner > 0 else 0.0
        v_outer_post  = (m_outer - m_vap) / (rho2 * A_outer) if A_outer > 0 else 0.0
        self.res['v_outer_annulus']   = round(v_outer,      2)
        self.res['v_inner_annulus']   = round(v_inner,      2)
        self.res['v_outer_post_vap']  = round(v_outer_post, 2)
        self.res['f_outer_feed_used'] = f_solved


def print_report(res, inputs):
    W = 65
    def header(title):
        print(f"\n{'─'*W}")
        print(f"  {title}")
        print(f"{'─'*W}")
    def row(label, value, unit="", note=""):
        line = f"  {label:<32} {value:>14}  {unit}"
        if note:
            line += f"   [{note}]"
        print(line)
    def flag(condition, msg):
        if condition:
            print(f"  ⚠  {msg}")

    print()
    print("=" * W)
    print("   REVERSE-FLOW ANNULAR MICRO-JET COMBUSTOR DESIGN TOOL")
    print("   KJ66-Style Architecture  |  Rev 21 — Dome Flux, Chamber Length, Stability Fixes")
    print("=" * W)

    oc_mm  = res['casing_od_mm']
    ol_mm  = res['outer_liner_od_mm']
    il_mm  = res['inner_liner_od_mm']
    st_mm  = res['shaft_tunnel_od_mm']
    scale  = oc_mm / 30.0

    def bar(od, id_=None):
        half_chars = int(round(oc_mm / 2 / scale))
        r_od = int(round((od / 2) / scale))
        r_id = int(round((id_ / 2) / scale)) if id_ is not None else 0
        line = " " * (half_chars - r_od) + "│" + "█" * (r_od - r_id) + " " * (r_id * 2) + "█" * (r_od - r_id) + "│"
        return line

    print()
    print("  RADIAL CROSS-SECTION (schematic):")
    print(f"  OD={oc_mm:.1f}mm  OL={ol_mm:.1f}mm  IL={il_mm:.1f}mm  ST={st_mm:.1f}mm")
    print()
    print("  " + "·" * 32 + " ← OUTER CASING")
    print("  " + bar(oc_mm, ol_mm) + " ← OUTER ANNULUS (feed air)")
    print("  " + bar(ol_mm, il_mm) + " ← COMBUSTION ANNULUS")
    print("  " + bar(il_mm, st_mm) + " ← INNER ANNULUS  (feed air)")
    print("  " + bar(st_mm) +        " ← SHAFT TUNNEL")
    print()

    header("[1] THERMODYNAMICS")
    row("Inlet Total Pressure P2", f"{res['P2_Pa']/1000:.1f}", "kPa")
    row("Inlet Total Temperature T2", f"{res['T2_K']:.0f}", "K")
    row("Target Turbine Inlet Temp", f"{inputs['target_tit_k']:.0f}", "K")
    clamp_note = "CLAMPED — real value higher" if res['T_primary_zone_clamped'] else "rich φ=1.6"
    row("Est. Primary Zone Temperature", f"{res['T_primary_zone_est']:.0f}", "K", clamp_note)
    if res['T_primary_zone_clamped']:
        print(f"  (Raw calculated: {res['T_primary_zone_raw']:.0f} K — clamped to 1950 K  [φ=1.6 rich zone cap])")
    row("Avg Cp Used", f"{res['cp_used']:.1f}", "J/kg·K")
    row("Inlet Air Density ρ2", f"{res['rho2']:.3f}", "kg/m³")

    header("[2] FLOW BUDGET")
    m_total = res['mdot_air']
    row("Total Air Mass Flow", f"{m_total:.4f}", "kg/s")
    row("Total Fuel Mass Flow", f"{res['mdot_fuel']*1000:.2f}", "g/s")
    row("Overall AFR", f"{res['overall_AFR']:.1f}", "kg_air/kg_fuel")
    row("Overall Equivalence Ratio φ", f"{res['overall_phi']:.3f}", "")
    hlf = res.get('heat_loss_factor', 1.00)
    print()
    if abs(hlf - 1.00) > 0.005:
        print(f"  NOTE: heat_loss_factor={hlf:.2f} applied — fuel flow increased by {(1/hlf-1)*100:.0f}%")
        print(f"  (More fuel needed to compensate for casing heat loss)")
    else:
        print(f"  NOTE: Adiabatic model (heat_loss_factor=1.00) — no casing loss modelled.")
        print(f"  KJ66 design-point phi=0.28 is correct at full-thrust TIT=1123K.")
        print(f"  Schreckling 3.8 g/s ref is part-throttle (~1040K TIT), not a sizing target.")
    print()
    print("  Air Split:")
    row("  Film Cooling (Effusion)", f"{res['m_film_cooling']/m_total:.1%}", "", f"{res['m_film_cooling']*1000:.1f} g/s  ({res['film_n_rows']} rows)")
    row("  Primary — Vaporizer Scoops", f"{res.get('split_primary_vap',0)/m_total:.1%}", "", f"{res.get('split_primary_vap',0)*1000:.1f} g/s")
    row("  Primary — Liner Holes", f"{res['split_primary_liner']/m_total:.1%}", "", f"{res['split_primary_liner']*1000:.1f} g/s")
    row("  Secondary Zone", f"{res['split_secondary']/m_total:.1%}", "", f"{res['split_secondary']*1000:.1f} g/s  φ→0.60")
    row("  Dilution Zone", f"{res['split_dilution']/m_total:.1%}", "", f"{res['split_dilution']*1000:.1f} g/s")
    phi_pri = res.get('phi_primary_actual', 1.6)
    row("  Primary zone φ_primary", f"{phi_pri:.2f}", "",
        "all primary air (liner + scoops)  [target: 1.0–1.6]")
    flag(phi_pri < 0.9,
         f"φ_primary = {phi_pri:.2f} < 0.9 — primary zone too lean. "
         "Blowout risk. Reduce primary air fraction or increase fuel.")
    flag(phi_pri > 2.0,
         f"φ_primary = {phi_pri:.2f} > 2.0 — primary zone very rich. "
         "Excessive soot/carbon formation. Increase primary air.")
    flag(res['split_dilution'] <= 0, "ZERO dilution air — check overall equivalence ratio.")

    header("[3] FOUR-WALL ANNULAR GEOMETRY")
    print("  OUTER BOUNDARY:")
    row("  Outer Casing OD", f"{res['casing_od_mm']:.2f}", "mm")
    row("  Outer Casing ID", f"{res['casing_id_mm']:.2f}", "mm")
    row("  Outer Feed Annulus Gap", f"{res['outer_annulus_gap_mm']:.2f}", "mm")
    row("  Outer Annulus Velocity", f"{res['v_outer_annulus']:.1f}", "m/s")
    print()
    print("  OUTER LINER:")
    print(f"  Liner Material: {res['liner_material_name']}  (α = {res['liner_alpha_per_K']*1e6:.1f}e-6 /K)")
    row("  Outer Liner OD (HOT / nominal)", f"{res['outer_liner_od_hot_mm']:.3f}", "mm", "reference only — do NOT use in CAD")
    row("  Outer Liner OD (COLD BUILD)",    f"{res['outer_liner_od_cold_mm']:.3f}", "mm", "← use this dimension in CAD / NX")
    row("  Thermal Radial Offset",          f"{res['outer_liner_thermal_offset_mm']:.3f}", "mm", "subtracted for cold-build fit")
    flag(res.get('thermal_offset_warn_outer', False),
         f"Outer liner thermal offset = {res['outer_liner_thermal_offset_mm']:.3f} mm "
         f"(raw: {res.get('outer_liner_thermal_offset_unclamped_mm', 0):.3f} mm). "
         "Large expansion — verify liner material rating and TIT. "
         "Cold-build dimension subtracts full computed offset.")
    row("  Outer Liner ID", f"{res['outer_liner_id_mm']:.2f}", "mm")
    print()
    print("  COMBUSTION ANNULUS:")
    row("  Radial Height H_comb", f"{res['combustion_gap_mm']:.2f}", "mm")
    row("  Mean Diameter D_mean", f"{res['D_mean_comb_mm']:.2f}", "mm")
    row("  Cross-Section Area", f"{res['combustion_annulus_A']*1e6:.1f}", "mm²")
    row("  Total Length", f"{res['chamber_length_mm']:.1f}", "mm")
    row("  L / D_mean ratio", f"{res['chamber_length_mm']/res['D_mean_comb_mm']:.2f}", "")
    print()
    print("  Zone Lengths:")
    row("  Primary Zone", f"{res['L_primary_mm']:.1f}", "mm")
    row("  Secondary Zone", f"{res['L_secondary_mm']:.1f}", "mm")
    row("  Dilution Zone", f"{res['L_dilution_mm']:.1f}", "mm")
    print()
    print("  INNER LINER:")
    row("  Inner Liner OD (HOT / nominal)", f"{res['inner_liner_od_hot_mm']:.3f}", "mm", "reference only — do NOT use in CAD")
    row("  Inner Liner OD (COLD BUILD)",    f"{res['inner_liner_od_cold_mm']:.3f}", "mm", "← use this dimension in CAD / NX")
    row("  Thermal Radial Offset",          f"{res['inner_liner_thermal_offset_mm']:.3f}", "mm", "subtracted for cold-build fit")
    flag(res.get('thermal_offset_warn_inner', False),
         f"Inner liner thermal offset = {res['inner_liner_thermal_offset_mm']:.3f} mm "
         f"(raw: {res.get('inner_liner_thermal_offset_unclamped_mm', 0):.3f} mm). "
         "Significant expansion — verify diametral clearance to shaft tunnel.")
    row("  Inner Liner ID", f"{res['inner_liner_id_mm']:.2f}", "mm")
    print()
    print("  INNER BOUNDARY:")
    row("  Shaft Tunnel OD", f"{res['shaft_tunnel_od_mm']:.2f}", "mm")
    row("  Shaft Tunnel ID", f"{res['shaft_tunnel_id_mm']:.2f}", "mm")
    row("  Inner Feed Annulus Gap", f"{res['inner_annulus_gap_mm']:.2f}", "mm")
    row("  Inner Annulus Velocity", f"{res['v_inner_annulus']:.1f}", "m/s")
    flag(res['v_outer_annulus'] > 50, f"Outer annulus velocity {res['v_outer_annulus']:.1f} m/s is high.")
    flag(res['v_inner_annulus'] > 50, f"Inner annulus velocity {res['v_inner_annulus']:.1f} m/s is high.")
    flag(res['combustion_gap_mm'] < 8, f"Combustion annulus height {res['combustion_gap_mm']:.1f} mm is very narrow.")

    f_split    = res.get('f_outer_solved', res.get('f_outer_feed_used', 0.65))
    converged  = res.get('split_converged', False)
    conv_label = "pressure-balance solved" if converged else "⚠ FALLBACK — minimised |ΔP| (no zero-crossing)"
    print(f"\n  Feed split: {f_split:.1%} outer / {1-f_split:.1%} inner  [{conv_label}]")
    if converged:
        row("  ΔP outer path", f"{res['dP_outer_path_Pa']:.1f}", "Pa")
        row("  ΔP inner path", f"{res['dP_inner_path_Pa']:.1f}", "Pa")
        row("  Solver tolerance", f"{res['split_solver_tol_Pa']:.1f}", "Pa", "= 1% of mean path ΔP")
        row("  Iterations", f"{res.get('split_solver_iters', '—')}", "",
            f"residual = {res.get('split_solver_residual_Pa', '—')} Pa")
    else:
        print("  ⚠ Pressure paths cannot be balanced — annulus velocities are mismatched.")
        print("  Golden-section minimiser found the closest-to-balanced split.")
        print(f"  Final residual: {res.get('split_solver_residual_Pa', '—')} Pa  "
              f"after {res.get('split_solver_iters', '—')} iterations.")
        print("  Consider matching outer/inner annulus velocities for a true pressure balance.")
        if res.get('split_solver_error'):
            print(f"  Solver exception: {res['split_solver_error']}")

    header("[4] VAPORIZER TUBES  (\"Candy Canes\"  — Air-Assist)")
    row("Number of Tubes (snapped)", f"{res['vap_n']}", "", f"raw={res['vap_n_raw']:.1f}, from {[6,8,12]}")
    row("Scoop Inlet Diameter", f"{res['vap_scoop_dia_mm']:.2f}", "mm", f"@ {res['vap_scoop_vel_actual']:.1f} m/s  (bell-mouth, 15–25° flare)")
    row("Tube OD / Bore ID", f"{res['vap_od_mm']:.2f} / {res['vap_id_mm']:.2f}", "mm")
    row("Crimp Orifice ID (fabrication)", f"{res['vap_crimp_dia_mm']:.2f}", "mm", "← drill/ream at 180° bend tip  [ref: 2–3.5 mm, S5]")
    row("U-bend Radius (1.5×OD)", f"{res['vap_bend_radius_mm']:.1f}", "mm")
    print()
    print("  Per-tube mass flow (air-assist):")
    row("  Fuel per tube", f"{res['vap_m_fuel_per_tube']:.2f}", "g/s")
    row("  Air per tube (scoop)", f"{res['vap_m_air_per_tube']:.2f}", "g/s")
    row("  Tube AFR", f"{res['vap_AFR_tube']:.2f}", "", "φ_tube from primary air fraction")
    print()
    print("  Outer annulus velocity accounting:")
    row("  Entry velocity (pre-scoop)", f"{res['v_outer_annulus']:.1f}", "m/s")
    row("  Post-scoop velocity", f"{res['v_outer_post_vap']:.1f}", "m/s", "after vaporizer extraction")
    print()
    print("  Fuel Pressure Drop (per tube):")
    row("  Tube gas flow regime", f"  {res['vap_tube_flow_regime']}", "",
        f"Re_gas = {res['fuel_Re_tube']:.0f}  (mix at 600K, {res['vap_mix_velocity']:.0f} m/s)")
    row("  Tube friction ΔP", f"{res['fuel_dP_tube_Pa']/1e5:.4f}", "bar")
    row("  Crimp nozzle ΔP", f"{res['fuel_dP_crimp_Pa']/1e5:.4f}", "bar")
    row("  TOTAL fuel ΔP", f"{res['fuel_dP_total_bar']:.3f}", "bar",
        f"reference target = {res['fuel_dP_target_bar']:.2f} bar  [S1: 0.4–0.8 bar]")
    flag(res['vap_id_mm'] < 3.5,
         f"Tube bore {res['vap_id_mm']:.1f} mm — consider reducing target_vap_mix_vel "
         "to widen bore and reduce coking risk.")
    flag(res['vap_id_mm'] > 7.0,
         f"Tube bore {res['vap_id_mm']:.1f} mm — verify OD fits within combustion annulus "
         "(OD/H_comb < 0.50 check below).")
    flag(res['fuel_dP_total_bar'] < 0.3, f"Fuel ΔP = {res['fuel_dP_total_bar']:.3f} bar is critically low.")
    flag(res['fuel_dP_total_bar'] > 1.2, f"Fuel ΔP = {res['fuel_dP_total_bar']:.3f} bar is high.")
    print()
    print("  Tube Geometry Validation:")
    L_tube_report = res.get('chamber_length_mm', 55) * 1.5 + 2 * res['vap_od_mm']
    od_over_H = res['vap_od_mm'] / res['combustion_gap_mm']
    LD_report  = L_tube_report / res['vap_id_mm']
    row("  Effective tube length", f"{L_tube_report:.0f}", "mm", "target: 60–120 mm")
    row("  Tube OD / H_comb", f"{od_over_H:.2f}", "",       "target: < 0.50 (tube fits in annulus)")
    row("  Tube L/D", f"{LD_report:.1f}", "",               "target: 8–30")
    flag(L_tube_report < 60, f"Tube length {L_tube_report:.0f} mm < 60 mm — may not complete vaporization.")
    flag(L_tube_report > 120, f"Tube length {L_tube_report:.0f} mm > 120 mm — cracking and thermal stress risk.")
    flag(od_over_H > 0.50, f"Tube OD/H = {od_over_H:.2f} — tube may not clear dome geometry.")
    flag(LD_report > 30, f"L/D = {LD_report:.1f} > 30 — excessive tube length risks thermal fatigue and coking.")
    print()
    print("  NOTE: Tube bore sized from gas-phase mixture velocity (85 m/s at 600K).")
    print("  High V_mix (80–100 m/s) creates turbulent mixing and atomisation needed")
    print("  for stable combustion. τ_vap ≈ 1 ms is sufficient for air-assist vaporizers.")

    header("[5] LINER HOLE SIZING  (Dual-Wall Annular)")
    print(f"  Hole area split — Outer: {res['hole_split_f_outer']:.1%}  Inner: {res['hole_split_f_inner']:.1%}")
    print(f"  (Proportional to liner circumference — Lefebvre method)")
    print()

    def hole_row(zone, qty_out, d_out, qty_in, d_in):
        print(f"  {zone:12s}  OUTER: {qty_out:2d} × {d_out:5.2f} mm     INNER: {qty_in:2d} × {d_in:5.2f} mm")

    hole_row("PRIMARY:",   res['pri_out_qty'], res['pri_out_mm'], res['pri_in_qty'],  res['pri_in_mm'])
    hole_row("SECONDARY:", res['sec_out_qty'], res['sec_out_mm'], res['sec_in_qty'],  res['sec_in_mm'])
    hole_row("DILUTION:",  res['dil_out_qty'], res['dil_out_mm'], res['dil_in_qty'],  res['dil_in_mm'])

    total_holes = (res['pri_out_qty'] + res['pri_in_qty'] + res['sec_out_qty'] + res['sec_in_qty'] + res['dil_out_qty'] + res['dil_in_qty'])
    print(f"\n  Total primary/secondary/dilution holes (both liners): {total_holes}")
    print()
    print("  Hole Mach Numbers (outer liner):")
    print(f"  Pri: Ma={res['Ma_hole_pri']:.3f}  Sec: Ma={res['Ma_hole_sec']:.3f}  Dil: Ma={res['Ma_hole_dil']:.3f}")
    flag(res['hole_compressible_warning'], "One or more hole zones have Ma > 0.3.")

    header("[5b] LINER FILM COOLING PATTERN  (Effusion Holes)")
    print(f"  {res['film_n_rows']} total rows: {res['film_n_rows_pri']} primary / {res['film_n_rows_sec']} secondary / {res['film_n_rows_dil']} dilution")
    row("Outer liner — holes per row", f"{res['film_holes_per_row_outer']}", "holes")
    row("Inner liner — holes per row", f"{res['film_holes_per_row_inner']}", "holes")
    row("Total effusion holes (both liners)", f"{res['film_total_holes']}", "")
    row("Total film area", f"{res['film_total_area_mm2']:.1f}", "mm²")

    header("[6] JET PENETRATION ANALYSIS  (Dilution Zone)")
    print(f"  Combustion annulus height H = {res['H_comb_mm']:.1f} mm")
    row("  Outer liner momentum flux J", f"{res['dil_J_outer']:.2f}", "")
    row("  Outer jet penetration depth", f"{res['dil_pen_outer_mm']:.1f}", "mm", f"{res['dil_pen_norm_outer']:.2f}×H  (ideal: 0.50×H)")
    row("  Inner liner momentum flux J", f"{res['dil_J_inner']:.2f}", "")
    row("  Inner jet penetration depth", f"{res['dil_pen_inner_mm']:.1f}", "mm", f"{res['dil_pen_norm_inner']:.2f}×H  (ideal: 0.50×H)")
    flag(res['dil_pen_norm_outer'] < 0.35, f"Outer dil. jets only penetrate {res['dil_pen_norm_outer']:.2f}×H.")
    flag(res['dil_pen_norm_outer'] > 0.70, f"Outer dil. jets penetrate {res['dil_pen_norm_outer']:.2f}×H. Combined coverage = {res['dil_combined_cov']:.2f}×H.")
    flag(res['dil_pen_norm_inner'] < 0.35, f"Inner dil. jets only penetrate {res['dil_pen_norm_inner']:.2f}×H.")

    header("[7] LINER STRUCTURAL CHECKS")
    row("  Outer liner hoop stress", f"{res['sigma_hoop_outer_MPa']:.1f}", "MPa")
    row("  Inner liner hoop stress", f"{res['sigma_hoop_inner_MPa']:.1f}", "MPa")
    row("  Allowable — 304 SS @ 900°C", f"{res['sigma_allow_ss304_MPa']:.0f}", "MPa")
    row("  Allowable — Inconel 625 @ 900°C", f"{res['sigma_allow_in625_MPa']:.0f}", "MPa")
    ok_ss = "✓ OK" if res['liner_ss304_ok'] else "✗ EXCEEDS"
    ok_in = "✓ OK" if res['liner_in625_ok'] else "✗ EXCEEDS"
    print(f"  304 SS:       {ok_ss}  |  Inconel 625:  {ok_in}")
    flag(not res['liner_in625_ok'], "Hoop stress exceeds Inconel 625 allowable.")
    flag(not res['liner_ss304_ok'] and res['liner_in625_ok'], "Hoop stress exceeds 304 SS allowable but within Inconel 625 limits.")
    print()
    row("  Outer liner ligament", f"{res['ligament_outer_mm']:.1f}", "mm")
    row("  Inner liner ligament", f"{res['ligament_inner_mm']:.1f}", "mm")
    row("  Required minimum", f"{res['ligament_min_req_mm']:.1f}", "mm")
    flag(not res['ligament_outer_ok'], f"Outer liner dilution holes too close — {res['ligament_outer_mm']:.1f} mm ligament.")
    flag(not res['ligament_inner_ok'], f"Inner liner dilution holes too close — {res['ligament_inner_mm']:.1f} mm ligament.")
    if res['ligament_outer_ok'] and res['ligament_inner_ok']:
        print("  ✓ Ligament spacing acceptable on both liners.")

    header("[8] COMBUSTION LOADING PARAMETER  (Lefebvre stability criterion)")
    row("Inlet Pressure P3",           f"{res['CLP_P3_kPa']:.1f}", "kPa")
    row("Primary Zone Volume",         f"{res['CLP_V_primary_m3']*1e6:.1f}", "cm³")
    row("Total Chamber Volume",        f"{res['CLP_V_total_m3']*1e6:.1f}", "cm³")
    row("CLP  (simple)",               f"{res['CLP']:.2f}", "kg/(s·kPa·m³)")
    row("CLP_T (temp-corrected)",      f"{res['CLP_T']:.2f}", "kg/(s·kPa·m³)")
    row("Combustor residence time τ",  f"{res['tau_comb_ms']:.2f}", "ms",
        "at mean comb. temp — target: >2 ms  |  flag <1.5 ms")
    row("Combustion efficiency η",     f"{res['eta_comb_estimated']:.3f}", "",
        "CLP-based estimate  [Lefebvre Table 5.1]")
    if res['CLP_stable']:
        print("  ✓ CLP in stable range — good blowout margin.")
    elif res['CLP_acceptable']:
        print("  ⚠  CLP in acceptable range.")
    else:
        flag(True, f"CLP = {res['CLP']:.1f} EXCEEDS blowout limit (15).")
    flag(res['tau_comb_ms'] < 1.5,
         f"τ_comb = {res['tau_comb_ms']:.2f} ms < 1.5 ms (mean-temp basis) — blowout risk "
         "at reduced throttle. Increase chamber volume (longer liner or wider combustion gap).")

    header("[9] MIXING INDEX  (first-order, NOT a true TTQ prediction)")
    row("Mixing Index (MI)", f"{res['mixing_index']:.3f}", "", "target < 0.25")
    row("Outer jet mix effectiveness", f"{res['mix_eff_outer']:.1%}", "")
    row("Inner jet mix effectiveness", f"{res['mix_eff_inner']:.1%}", "")
    cov = res['dil_combined_cov']
    row("Dual-wall combined coverage", f"{cov:.2f}", "×H",
        "[<0.8: cold core  |  ~1.0: midplane met  |  1.0–2.0: good  |  >2.0: liner impingement]")
    if cov < 0.8:
        print("  ✗ Combined coverage < 0.8×H — jets don't reach mid-annulus.")
        print("    Hot core survives to turbine inlet. Increase hole size or count.")
    elif cov < 1.0:
        print("  ⚠  Combined coverage 0.8–1.0×H — jets approach but don't meet at midplane.")
        print("    Some radial temperature stratification expected.")
    elif cov <= 2.0:
        print(f"  ✓ Combined coverage {cov:.2f}×H — good radial mixing.")
        if cov < 1.2:
            print("    Jets just meeting at mid-annulus — ideal penetration balance.")
        elif cov <= 1.6:
            print("    Jets crossing mid-annulus slightly — good mixing, no impingement risk.")
        else:
            print("    Jets crossing well past mid-annulus — check opposite-wall ligament integrity.")
    else:
        print("  ⚠  Combined coverage > 2.0×H — each jet individually reaches the opposite liner.")
        print("    Risk of hot-spot damage at liner dilution holes. Increase dil_holes_per_vap.")
    if res['mixing_index'] > 0.25:
        print(f"  NOTE: MI = {res['mixing_index']:.3f} (single-wall Sturgess flag); combined coverage is primary metric.")

    print()
    print("=" * W)
    print("  DESIGN NOTES:")
    print("  • Liner wall = 0.5–1.5mm stainless steel (304 SS) or Inconel 625")
    print("  • All jet holes: deburred, edge-broken 0.1–0.2mm chamfer")
    print("  • Vaporizer tubes: 316 SS or Inconel 625 for thermal cycling resistance")
    print("  • Film cooling rows: first row immediately after primary holes; stagger 30°")
    print("  • Igniter plug: outer liner, primary zone, midway between vaporizers")
    print("=" * W)
    header("[10a] COMBUSTOR EXIT CONDITIONS  (turbine inlet)")
    
    row("Exit static pressure P4",      f"{res['exit_P4_kPa']:.2f}",   "kPa",
        f"P2 × (1 − {(1 - res['exit_P4_Pa']/res['P2_Pa'])*100:.0f}% ΔP)")
    row("Exit static temperature T4",   f"{res['exit_T4_K']:.0f}",     "K",
        "= target TIT (design point)")
    row("Exit total temperature T04",   f"{res['exit_T4_total_K']:.1f}","K")
    row("Exit total pressure P04",      f"{res['exit_P4_total_Pa']/1000:.2f}", "kPa")
    row("Exit gas density ρ4",          f"{res['exit_rho4']:.4f}",      "kg/m³")
    row("Exit mass flow (air + fuel)",  f"{res['exit_mdot_kg_s']:.4f}", "kg/s")
    row("Exit bulk velocity V4",        f"{res['exit_V4_m_s']:.1f}",    "m/s")
    row("Speed of sound a4",            f"{res['exit_a4_m_s']:.1f}",    "m/s",
        f"γ = {res['exit_gamma_hot']} (hot products)")
    row("Exit Mach number Ma4",         f"{res['exit_Ma4']:.4f}",       "")
    row("Cp at T4",                     f"{res['exit_cp4']:.1f}",       "J/kg·K")
    row("Specific enthalpy available",  f"{res['exit_h_avail_kJ_kg']:.1f}", "kJ/kg",
        "total enthalpy drop T04 → T_amb")
    flag(res['exit_Ma4_high'],
        f"Ma4 = {res['exit_Ma4']:.3f} > 0.25 — compressibility non-negligible. "
        "Isentropic stagnation values above are valid; static values are meaningfully "
        "lower than total. Feed total conditions (T04, P04) to turbine design.")
    flag(res['exit_Ma4'] > 0.50,
        f"Ma4 = {res['exit_Ma4']:.3f} — EXIT IS TRANSONIC-APPROACHING. "
        "Combustion annulus area may be too small for this mass flow and TIT.")
    print()
    print("  NOTE: P4 derived from design-target ΔP, not independently verified.")
    print(f"  γ = {res['exit_gamma_hot']} used here vs 1.4 used elsewhere in code.")
    print("  For turbine design use total conditions (T04, P04) — not static.")
    header("[10b] STABILITY CHECKS  (Recirculation · Flashback · Dome Heat)")
    print("  These checks flag fluid-dynamic failure modes invisible to static geometry sizing.")
    print("  A 1-D tool cannot guarantee stability — use as go/no-go screening only.")
    print()

    print("  CHECK 1 — Primary Jet Recirculation Strength:")
    row("  J_primary  (ρ_jet·Vjet² / ρ_hot·Vax²)", f"{res['stab_J_primary']:.1f}", "",
        "target: 5–80")
    row("  Primary jet velocity Vjet", f"{res['stab_V_jet_pri']:.1f}", "m/s")
    row("  Combustion zone axial Vax", f"{res['stab_V_axial']:.1f}", "m/s")
    row("  Hot-gas density ρ_hot", f"{res['stab_rho_pri']:.3f}", "kg/m³",
        f"at T_primary={res['T_primary_zone_est']:.0f} K")
    if res['stab_J_ok']:
        print("  ✓ J in stable range — recirculation vortex should sustain.")
    flag(res['stab_J_collapse'],
         f"J = {res['stab_J_primary']:.1f} < 5 — vortex collapse risk. "
         "Primary jets may be swept downstream: increase primary hole area or reduce annulus bulk velocity.")
    flag(res['stab_J_overpenetrate'],
         f"J = {res['stab_J_primary']:.1f} > 80 — jets likely over-penetrate and split the vortex. "
         "Reduce primary hole velocity: add more holes or enlarge diameter.")
    flag(res.get('stab_V_axial_high', False),
         f"Primary zone axial velocity V_axial = {res['stab_V_axial']:.1f} m/s > 30 m/s — "
         "recirculation vortex may collapse under high bulk flow. "
         "Increase combustion annulus cross-section or reduce annulus air velocity target.")

    print()
    print("  CHECK 2 — Vaporizer Flashback & Residence Time:")
    row("  Mixture velocity in tube", f"{res['stab_V_mix']:.1f}", "m/s",
        "threshold: >4.0 safe, <0.8 critical")
    row("  Tube L/D", f"{res['stab_LD_tube']:.1f}", "",
        "≥8 recommended for quench margin")
    row("  Vaporizer residence time τ", f"{res['stab_tau_vap_ms']:.2f}", "ms",
        "air-assist arch. threshold: >0.3 ms")
    if not res['stab_flashback_critical'] and not res['stab_flashback_risk'] and res['stab_LD_ok']:
        print("  ✓ Flashback margin adequate.")
    elif not res['stab_flashback_critical'] and not res['stab_flashback_risk']:
        print("  ✓ Mixture velocity adequate — but check L/D.")
    flag(res['stab_flashback_critical'],
         f"CRITICAL: V_mix = {res['stab_V_mix']:.1f} m/s — IMMEDIATE FLASHBACK RISK. "
         "Flame will propagate upstream into the vaporizer tube. Reduce tube bore or increase mass flow per tube.")
    flag(res['stab_flashback_risk'] and not res['stab_flashback_critical'],
         f"V_mix = {res['stab_V_mix']:.1f} m/s < 4.0 m/s — flashback possible under transients. "
         "Add a quench/choke orifice near the tube dome end.")
    flag(not res['stab_LD_ok'],
         f"Tube L/D = {res['stab_LD_tube']:.1f} < 8 with bore > 6mm — "
         "quench distance may be insufficient to extinguish a flame front inside the tube.")
    flag(res['stab_tau_vap_ms'] < 0.3,
         f"τ_vap = {res['stab_tau_vap_ms']:.2f} ms < 0.3 ms — tube may be too short for "
         "adequate fuel/air mixing before the 180° bend. Risk of coking or stratified combustion. "
         "Increase effective tube length or reduce mixture velocity.")

    print()
    print("  CHECK 3 — Dome Heat Flux Estimate (Convective + Radiative):")
    row("  Assumed dome wall temperature", f"{res['stab_T_wall_dome_K']:.0f}", "K",
        "= 80% TIT (unshielded dome conservative)")
    row("  Convective flux q_conv", f"{res['stab_q_conv_kW_m2']:.0f}", "kW/m²",
        f"Nu=0.037·Re^0.8·Pr^0.33  Re={res['stab_Re_imp']:.0f}  h={res['stab_h_conv']:.0f} W/m²K")
    row("  Radiative flux q_rad", f"{res['stab_q_rad_kW_m2']:.0f}", "kW/m²",
        f"ε={res['stab_epsilon_flame']:.2f}  F_view={res['stab_F_view']:.2f}  (grey-body, Lefebvre §9.3)")
    row("  TOTAL dome heat flux q", f"{res['stab_q_dome_kW_m2']:.0f}", "kW/m²",
        "[<300: OK  |  300–600: HIGH  |  >600: CRITICAL]")
    if res['stab_dome_ok']:
        print("  ✓ Dome heat flux acceptable.")
    flag(res['stab_dome_high'],
         f"q_dome = {res['stab_q_dome_kW_m2']:.0f} kW/m² — HIGH. "
         "Add a dome cooling air ring (≥5% total air) or transpiration cooling film.")
    flag(res['stab_dome_critical'],
         f"q_dome = {res['stab_q_dome_kW_m2']:.0f} kW/m² — CRITICAL. "
         "Dome will fail without active cooling. Add impingement-cooled dome insert "
         "or increase primary air fraction to lower T_flame.")
    print()
    print("  Methodology: convection via turbulent impingement (Dittus-Boelter variant).")
    print("  Radiation via grey-body model with phi-dependent emissivity and F_view=0.60.")
    print("  Treat as order-of-magnitude screening — CFD required for final clearance.")

    print()
    print("  CHECK 4 — Vaporizer Tube Thermal & Coking Risk:")
    print("  1-D resistance model: recirculation zone → tube wall → mixture")
    print(f"  T_primary (outside tube) = {res.get('T_primary_zone_est', 1950):.0f} K")
    print(f"  T_mixture (inside tube)  ≈ 600 K  (conservative estimate)")
    row("  Tube wall h_inner (Dittus-Boelter)", f"{res['stab_h_inner_vap']:.0f}", "W/m²K",
        "turbulent mixture flow in 4 mm bore")
    row("  T_wall outer face", f"{res['stab_T_vap_wall_outer_K']:.0f}", "K",
        f"({res['stab_T_vap_wall_outer_K']-273:.0f}°C)  — Inconel 625 cont. limit ~1100°C")
    row("  T_wall inner face", f"{res['stab_T_vap_wall_inner_K']:.0f}", "K",
        f"({res['stab_T_vap_wall_inner_K']-273:.0f}°C)  — Jet-A coking onset ~650 K (377°C)")
    if res['stab_vap_coking_risk']:
        print("  ⚠  COKING RISK: Inner wall > 650 K — Jet-A thermal cracking will deposit carbon.")
        print("     This is expected and IS the dominant life-limiting failure mode on real KJ66s.")
        print("     Mitigation: use 316SS or Inconel 625 tubes, inspect/replace every 10–20 hrs.")
        print("     Increasing air-assist flow (primary_air_vaporizer_fraction) reduces wall temp.")
    else:
        print("  ✓ Inner wall below Jet-A coking threshold.")
    if res['stab_vap_structural_risk']:
        print("  ⚠  STRUCTURAL RISK: Outer wall > 1100 K — 304SS will creep/fail rapidly.")
        print("     MATERIAL REQUIREMENT: Use Inconel 625 or 316SS for vaporizer tubes.")
    else:
        print("  ✓ Outer wall within Inconel 625 continuous service limit.")
    print()
    print("  NOTE: h_outer = 300 W/m²K assumed (recirculation zone impingement).")
    print("  Actual T_wall depends on local recirculation velocity and dome cooling.")

    header("[11] CAD GEOMETRY EXPORT  (copy-paste into CAD)")
    import json
    cad = res['cad_geometry']
    print(json.dumps(cad, indent=2))
    print()


user_inputs_6in = {
    'casing_od_inch':        6.00,
    'shaft_tunnel_od_inch':  1.65,
    'wall_thickness_mm':     1.50,
    'pressure_ratio':        1.50,
    'compressor_efficiency': 0.94,
    'mass_flow_air_kg_s':    0.6, # from .487
    'target_tit_k':          900.0,
    'liner_material':       '316SS',
}

kj66_inputs = {
    'casing_od_inch':        4.33,
    'shaft_tunnel_od_inch':  1.18,
    'wall_thickness_mm':     0.50,
    'pressure_ratio':        2.20,
    'compressor_efficiency': 0.74,
    'mass_flow_air_kg_s':    0.23,
    'target_tit_k':         1123.0,
    'liner_material':       '304SS',
    # KJ66-specific empirical calibration: inner annulus runs at 20 m/s
    'target_inner_annulus_vel': 20.0,
    # tau_min_s calibrated to reproduce KJ66 empirical chamber length ~65 mm.
    # With bulk-density V_ref (correct physics), default tau_min=0.002s gives 95 mm —
    # physically justified for a new design but 47% longer than the historical hardware.
    # 0.0014s with bulk density → ~67 mm, matching Schreckling [S1] chamber geometry.
    # General (non-KJ66) designs should keep the default tau_min_s=0.002s.
    'tau_min_s': 0.0014,
}

KEYS = [
    'casing_od_inch',
    'shaft_tunnel_od_inch',
    'wall_thickness_mm',
    'pressure_ratio',
    'compressor_efficiency',
    'mass_flow_air_kg_s',
    'target_tit_k',
]

if __name__ == "__main__":
    print()
    print("  REVERSE-FLOW ANNULAR COMBUSTOR DESIGN TOOL")
    print()
    print("  Presets available:")
    print("    1 — 6\" custom engine  (PR 1.5,  0.6 kg/s,  TIT 900 K)")
    print("    2 — KJ66 reference    (PR 2.2,  0.230 kg/s, TIT 1123 K)")
    print("    3 — Manual input")
    print()
    choice = input("  Select [1/2/3]: ").strip()

    if choice == "1":
        chosen_inputs = user_inputs_6in
        print("\n  Using: 6-inch engine preset")
    elif choice == "2":
        chosen_inputs = kj66_inputs
        print("\n  Using: KJ66 reference preset")
    else:
        print()
        print("  Enter values in this order (space-separated):")
        print("  casing_od_inch | shaft_tunnel_od_inch | wall_thickness_mm |")
        print("  pressure_ratio | compressor_efficiency | mass_flow_air_kg_s | target_tit_k")
        print()
        print("  Example (KJ66):  4.33 1.18 0.5 2.2 0.74 0.23 1123")
        while True:
            raw = input("  Enter: ").strip().split()
            if len(raw) != len(KEYS):
                print(f"  Expected {len(KEYS)} values, got {len(raw)}. Try again.")
                continue
            try:
                parsed = [float(v) if v.lower() != 'none' else None for v in raw]
                break
            except ValueError:
                print("  Invalid value — ensure all entries are numbers. Try again.")
        chosen_inputs = dict(zip(KEYS, parsed))
        mat_raw = input("  Material [304SS]: ").strip().upper() or '304SS'
        chosen_inputs['liner_material'] = mat_raw

    print(f"\n  Inputs: {chosen_inputs}")

    try:
        model   = MicroJetCombustor(chosen_inputs)
        results = model.run()
        print_report(results, chosen_inputs)
    except ValueError as e:
        print(f"\n  ✗ GEOMETRY ERROR: {e}")
        print("  Adjust casing OD, shaft tunnel OD, or wall thickness and retry.\n")

    print()
    while True:
        close = input("  Press Enter, X, or Q to exit: ").strip().lower()
        if close in ("", "x", "q", "5"):
            break



#  AIR FLOW PATH (Reverse-Flow):
#   Compressor exit → splits into OUTER and INNER annuli →
#   travels REARWARD alongside combustion zone →
#   enters combustion annulus via liner holes (primary, secondary, dilution) →
#   vaporizer tubes at rear dome inject fuel FORWARD (counter-flow) →
#   combustion products flow FORWARD toward turbine inlet
#
#  HOLE SPLIT METHODOLOGY (Lefebvre, Gas Turbine Combustion 3rd Ed.):
#   Total required hole area for each zone is split proportionally to
#   the circumference of the respective liner wall. Because the outer
#   liner has a larger circumference, it receives a greater share.
#   Split fraction = C_outer / (C_outer + C_inner) for outer liner
#   This ensures approximately equal jet momentum flux from both walls,
#   targeting jet penetration to the mid-height of the combustion annulus.
#
#  JET PENETRATION (Sturgess / Lefebvre correlation):
#   y_max / H = 1.15 * (J)^0.5 * (d_j / H)
#   where J = momentum flux ratio = (rho_j * V_j^2) / (rho_g * V_g^2)
#   Target: y_max / H ≈ 0.5  (jet reaches combustion annulus midplane)