import math
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

class MicroJetCombustor:
    def __init__(self, inputs):
        self.inputs = inputs
        self.res = {}
        self.R = 287.05     # J/(kg·K) — specific gas constant for air
        self.GAMMA = 1.4    # ratio of specific heats for air

        self.DESIGN_PARAMS = {
            # --- Annulus Velocities ---
            'target_annulus_vel': 35.0,     # m/s  — air in feed annuli; low for uniform pressure
            'target_tube_liq_vel': 3.0,     # m/s  — liquid fuel in vaporizer tube bore
            'vap_pitch_mm': 35.0,           # mm   — circumferential pitch between vaporizer tubes
            # --- Combustor Sizing ---
            'target_pressure_drop': 0.04,   # 4%   — dP across liner (industry standard for micro-GT)
            'discharge_coeff_hole': 0.60,   # Cd   — sharp-edged holes (Lefebvre recommendation)
            # --- Zone Length Fractions of Total Liner Length ---
            'frac_primary_zone': 0.30,      # 30%  — primary/rich zone (flame stabilization)
            'frac_secondary_zone': 0.30,    # 30%  — secondary zone (soot burnout, CO oxidation)
            'frac_dilution_zone': 0.40,     # 40%  — dilution zone (temperature profile shaping)
            # --- Feed Annulus Air Split ---
            # This controls how total air mass flow is divided between the outer
            # feed annulus (casing-to-outer-liner gap) and the inner feed annulus
            # (inner-liner-to-shaft-tunnel gap).
            #
            # Critical:
            # In real hardware, the mass split is determined by equalization of
            # static pressure at the liner entry, which depends on hydraulic
            # diameter, friction length, turn losses, and vaporizer extraction
            # from the outer annulus. A full solution requires iterating on
            # pressure loss in each path.
            #
            # This tool instead ASSUMES mass split = area split (f_outer_feed
            # applied to both geometry and mass flow simultaneously), which
            # produces equal bulk velocities in both annuli by construction.
            # That is a valid first-order sizing assumption for a conceptual
            # design tool, but is not physically derived.
            #
            # Typical range: 0.55–0.65 (outer annulus larger → more flow).
            # Outer annulus draws more because:
            #   (a) larger mean radius → larger area for same gap width
            #   (b) vaporizer tubes extract primary air from outer annulus
            # Adjust this parameter if your geometry shows velocity imbalance.
            'f_outer_feed': 0.60,
            # --- Hole Split Between Inner and Outer Liner ---
            # Computed dynamically from liner circumferences; this is an override if set > 0
            # Set to None to use circumference-proportional split (recommended)
            'outer_liner_hole_fraction': None,  # e.g. 0.60; None = auto from circumference ratio
            # --- Air Budget Fractions ---
            'film_cooling_fraction': 0.12,  # 12%  — reserved for liner film cooling / effusion
            # --- Vaporizer Air Ingestion ---
            # Fraction of total primary-zone air that enters via the vaporizer tubes
            # rather than through the primary liner holes.
            #
            # In a reverse-flow air-assist vaporizer (candy-cane style), the tube
            # dome scoop draws both liquid fuel AND a portion of primary air from the
            # outer feed annulus. This air mixes with fuel vapor inside the tube and
            # co-injects forward into the recirculation zone.
            #
            # Ignoring this causes three compounding errors:
            #   1. Primary liner hole area oversized (air already delivered via tubes)
            #   2. Outer annulus velocity underestimated (extraction not accounted for)
            #   3. Vaporizer exit velocity underpredicted (more mass, same exit area)
            #
            # Typical range: 0.15–0.35. Default 0.25 is mid-range for KJ66-class hardware.
            # Set to 0.0 to revert to fuel-only tube model (not physically correct).
            'primary_air_vaporizer_fraction': 0.25,
            # --- Vaporizer Tube Count Constraint ---
            # Raw count is computed from circumferential pitch at mean annulus diameter,
            # then snapped to the nearest value in this list.
            #
            # 6, 8, and 12 give symmetric angular patterns (60°, 45°, 30° spacing).
            # Symmetry matters for ignition propagation, thermal balance, and strut alignment.
            # Counts outside this list are snapped to nearest allowed value.
            'allowed_vaporizer_counts': [6, 8, 12],
            # --- Vaporizer Scoop Inlet (HIGH-SEVERITY — mandatory for correct air metering) ---
            # The scoop is the opening in the dome that draws primary air into the
            # vaporizer tube. Without a defined scoop geometry, m_air_vap is just a
            # number in the model with no hardware to deliver it.
            # Undersized scoop → over-rich primary bubble → coking → tube blockage
            # Oversized scoop → recirculation collapses → lean blowout / hard light-off
            'vap_scoop_target_vel_m_s': 22.0,  # m/s — low velocity for stable entrainment
            'vap_scoop_cd':             0.72,  # Cd — flared bell-mouth or open tube
            # --- Explicit Film-Cooling Pattern (HIGH-SEVERITY — mandatory for liner survival) ---
            # Without hole layout, "12% film air" exists only on paper.
            # Liner burnout typically occurs because film rows are too sparse or
            # placed too far downstream from primary holes.
            'film_row_pitch_mm':    38.0,  # mm  — axial spacing between cooling rows
            'film_hole_dia_mm':      1.20,  # mm  — discrete effusion hole diameter
            'film_rows_primary':        1,  # rows in primary zone
            'film_rows_secondary':      1,  # rows in secondary zone
            'film_rows_dilution':       2,  # rows in dilution zone (longest, most exposed)
        }

        # --- Fuel Properties (Jet-A / Kerosene) ---
        self.FUEL = {
            "LHV":       43.0e6,    # J/kg  — lower heating value
            "STOICH_AFR": 14.7,     # kg_air/kg_fuel — stoichiometric air-fuel ratio
            "RHO_LIQ":   800.0,     # kg/m³ — liquid density at ambient
            "T_BOIL":    450.0,     # K     — approximate vaporization temperature
        }

    # -------------------------------------------------------------------------
    def _get_cp(self, T_kelvin):
        """
        Variable specific heat for air/combustion mix. Linear fit vs T.

        Calibrated against NACA SP-273 / Gordon-McBride equilibrium data for
        lean-to-stoichiometric Jet-A combustion products at 1–4 atm:
          ~1005 J/kgK at 300 K rising to ~1150 J/kgK around 1500 K.

        Previous model (950 + 0.21T) overestimated Cp by 5–8% above 800 K,
        inflating fuel flow calculations. Corrected to 1005 + 0.10*(T-300),
        clamped to 1150 J/kgK max — physically consistent with lean combustion
        products at micro-jet pressure ratios (1.5–4 atm).
        """
        T  = max(300, min(T_kelvin, 2000))
        cp = 1005.0 + 0.10 * (T - 300)
        return min(cp, 1150.0)  # J/(kg·K) — hard cap at physically realistic maximum

    def thermodynamics(self):
        """Compressor exit conditions — air entering combustor."""
        P_amb = 101325   # Pa
        T_amb = 288.15   # K
        PR    = self.inputs['pressure_ratio']
        eta_c = self.inputs['compressor_efficiency']

        P2     = P_amb * PR
        T2_iso = T_amb * (PR ** ((self.GAMMA - 1) / self.GAMMA))
        T2     = T_amb + (T2_iso - T_amb) / eta_c   # isentropic + efficiency
        rho2   = P2 / (self.R * T2)

        self.res['P2_Pa']  = P2
        self.res['T2_K']   = T2
        self.res['rho2']   = rho2

    def mass_flow_and_fuel(self):
        """Determines air mass flow and required fuel flow for target TIT."""
        m_air = self.inputs.get('mass_flow_air_kg_s')
        if m_air is None:
            # Auto-scale from casing OD (reference: 6" case → 0.45 kg/s)
            ref_od   = 0.1524
            ref_flow = 0.45
            od_m     = self.inputs['casing_od_inch'] * 0.0254
            m_air    = ref_flow * (od_m ** 2 / ref_od ** 2)

        target_tit = self.inputs['target_tit_k']
        T2         = self.res['T2_K']
        T_avg      = (target_tit + T2) / 2
        cp_avg     = self._get_cp(T_avg)
        comb_eff   = 0.96

        energy_req = m_air * cp_avg * (target_tit - T2)
        m_fuel     = energy_req / (self.FUEL['LHV'] * comb_eff)

        self.res['mdot_air']    = m_air
        self.res['mdot_fuel']   = m_fuel
        self.res['overall_AFR'] = m_air / m_fuel
        self.res['overall_phi'] = self.FUEL['STOICH_AFR'] / self.res['overall_AFR']
        self.res['cp_used']     = cp_avg

    def zonal_analysis(self):
        """
        Splits combustion air between primary, secondary, and dilution zones.
        Also reserves a fraction for film/effusion cooling of the liner walls.

        Primary zone  φ ≈ 1.6  (rich, stable, avoids soot at too-lean conditions)
        Secondary zone target: drive overall φ down to ~0.6 (CO burnout region)
        Dilution zone: remainder (profile shaping for turbine)

        Film cooling air does NOT pass through liner holes — it enters through
        closely-spaced effusion holes or slot-cooling strips along the liner.
        It is accounted for FIRST before computing combustion-zone splits.
        """
        m_fuel  = self.res['mdot_fuel']
        stoich  = self.FUEL['STOICH_AFR']
        m_total = self.res['mdot_air']

        # Film cooling reservation
        f_film  = self.DESIGN_PARAMS['film_cooling_fraction']
        m_film  = m_total * f_film
        m_comb  = m_total - m_film   # air available for combustion zones

        # Primary zone air — total required at φ=1.6
        phi_pri       = 1.6
        m_air_pri_total = (m_fuel * stoich) / phi_pri

        # Feasibility check: primary zone cannot demand more air than is available
        if m_air_pri_total > m_comb:
            raise ValueError(
                f"Primary zone requires {m_air_pri_total*1000:.1f} g/s air but only "
                f"{m_comb*1000:.1f} g/s is available after film cooling reservation. "
                f"Reduce film_cooling_fraction, increase mass_flow_air, or lower target_tit_k."
            )

        # Split primary air: fraction delivered via vaporizer dome scoops vs liner holes.
        # The vaporizer fraction does NOT flow through primary liner holes — it is
        # extracted from the outer feed annulus at the dome and co-injected with fuel vapor.
        f_vap         = self.DESIGN_PARAMS['primary_air_vaporizer_fraction']
        m_air_vap     = m_air_pri_total * f_vap      # air through vaporizer tubes
        m_air_pri_liner = m_air_pri_total * (1.0 - f_vap)  # air through primary liner holes

        # m_air_pri (total) kept for equivalence ratio accounting
        m_air_pri = m_air_pri_total

        # Secondary zone — brings cumulative φ to ~0.6
        phi_sec_target = 0.6
        m_cumulative   = (m_fuel * stoich) / phi_sec_target
        m_air_sec      = m_cumulative - m_air_pri

        # Dilution zone — remainder of combustion air
        m_air_dil = m_comb - m_air_pri - m_air_sec
        if m_air_dil < 0:
            m_air_dil = 0
            m_air_sec = m_comb - m_air_pri

        self.res['m_film_cooling']      = m_film
        self.res['split_primary']       = m_air_pri          # total primary air (φ accounting)
        self.res['split_primary_liner'] = m_air_pri_liner    # primary air via liner holes only
        self.res['split_primary_vap']   = m_air_vap          # primary air via vaporizer scoops
        self.res['split_secondary']     = m_air_sec
        self.res['split_dilution']      = m_air_dil
        self.res['m_combustion_air']    = m_comb

        # Primary zone temperature estimate (first zone, rich combustion)
        # Uses energy balance over primary zone only (fuel burned at phi=1.6).
        # Physical upper bound ~2050 K: above this, mixing lag, quenching from film
        # cooling air ingestion, and finite flame speed all suppress actual temperature.
        cp_flame        = self._get_cp(1900)       # representative rich-zone Cp
        m_fuel_burned   = m_air_pri / stoich        # fuel that reacts in primary zone
        Q_zone1         = m_fuel_burned * self.FUEL['LHV']
        m_zone1_total   = m_air_pri + m_fuel        # all fuel enters + primary air
        T_pri_raw       = self.res['T2_K'] + Q_zone1 / (m_zone1_total * cp_flame)
        T_pri_clamped   = min(T_pri_raw, 2050.0)
        self.res['T_primary_zone_est']     = T_pri_clamped
        self.res['T_primary_zone_raw']     = T_pri_raw
        self.res['T_primary_zone_clamped'] = T_pri_raw > 2050.0   # flag if clamp active

    def mechanical_geometry(self):
        """
        Sizes the four-wall annular combustor geometry.

        BOUNDARIES (radial, outward → inward):
          1. Outer Casing     — user input (fixed)
          2. Outer Liner      — computed from outer annulus velocity target
          3. Inner Liner      — computed from annulus height and shaft tunnel
          4. Shaft Tunnel     — user input (fixed)

        DESIGN LOGIC:
          a) Total annular flow area required (both feed annuli) from velocity target.
          b) The outer feed annulus area = A_total - inner feed annulus area.
             Inner annulus constrained by shaft tunnel OD + minimum clearance.
          c) Combustion annulus area fills the space between the two liners.
          d) L/D uses mean combustion annulus diameter as characteristic dimension.
        """
        # --- Parse physical boundaries ---
        wall_m     = self.inputs['wall_thickness_mm'] / 1000.0
        od_casing  = self.inputs['casing_od_inch'] * 0.0254          # outer casing OD (m)
        id_casing  = od_casing - 2 * wall_m                          # outer casing ID (m)
        od_tunnel  = self.inputs['shaft_tunnel_od_inch'] * 0.0254    # shaft tunnel OD (m)
        id_tunnel  = od_tunnel - 2 * wall_m                          # shaft tunnel ID (m; structural)

        # --- Total annulus area needed for air feed at target velocity ---
        v_ann      = self.DESIGN_PARAMS['target_annulus_vel']
        rho2       = self.res['rho2']
        m_air      = self.res['mdot_air']
        A_feed_total = m_air / (rho2 * v_ann)   # combined outer + inner feed annuli area

        # --- Outer liner OD ---
        # Outer feed annulus: between casing ID and outer liner OD
        # Inner feed annulus: between inner liner ID and shaft tunnel OD
        # Split feed area proportionally; start by assuming equal split, then iterate.
        # For the KJ66-class geometry, the outer annulus is larger.
        # We use a first-order estimate: outer annulus gets (R_casing² - R_guess²) of area.

        # The total available cross-section for casing bore
        A_casing_id = math.pi * (id_casing / 2) ** 2
        A_tunnel_od = math.pi * (od_tunnel / 2) ** 2

        # --- Feed annulus area split ---
        # See DESIGN_PARAMS['f_outer_feed'] for full documentation of this assumption.
        # Short version: mass split is IMPOSED equal to area split — velocities will
        # be equal by construction, but this is an assumption, not a physics solution.
        f_outer_feed = self.DESIGN_PARAMS['f_outer_feed']
        A_outer_feed   = A_feed_total * f_outer_feed
        A_inner_feed   = A_feed_total * (1.0 - f_outer_feed)

        # Outer liner OD: casing ID area minus outer feed area
        r_casing_id    = id_casing / 2
        r_outer_liner_od = math.sqrt(r_casing_id ** 2 - A_outer_feed / math.pi)
        outer_liner_od = r_outer_liner_od * 2

        # Inner liner ID: shaft tunnel OD + inner feed annulus area
        r_tunnel_od    = od_tunnel / 2
        r_inner_liner_id = math.sqrt(r_tunnel_od ** 2 + A_inner_feed / math.pi)
        inner_liner_id = r_inner_liner_id * 2

        # Apply wall thickness to get the other liner surface
        outer_liner_id = outer_liner_od - 2 * wall_m
        inner_liner_od = inner_liner_id + 2 * wall_m

        # --- Sanity check: inner liner must be smaller than outer liner ---
        combustion_gap = (outer_liner_id - inner_liner_od) / 2   # radial height (m)
        if combustion_gap <= 0.005:  # must be at least 5mm clear
            raise ValueError(
                f"Geometry infeasible: combustion annulus gap = {combustion_gap*1000:.1f} mm. "
                f"Reduce shaft tunnel OD or increase casing OD."
            )

        # --- Combustion annulus area ---
        A_comb = math.pi * ((outer_liner_id / 2) ** 2 - (inner_liner_od / 2) ** 2)

        # --- Mean diameter of combustion annulus ---
        D_mean = (outer_liner_id + inner_liner_od) / 2   # characteristic diameter for L/D

        # --- Combustor length ---
        # Use L/D based on MEAN combustion annulus diameter.
        # For reverse-flow annular combustors, L/D_mean ≈ 1.5–2.0 is typical.
        max_LD = 1.65
        length = D_mean * max_LD

        # --- Verify actual annulus velocities ---
        # The outer annulus loses m_air_vap to the vaporizer dome scoops partway along
        # its length. Velocity reported here is the ENTRY velocity (before extraction),
        # which is the design-controlling value for pressure uniformity at the liner face.
        # Post-extraction velocity in the downstream half of the outer annulus will be
        # slightly lower — conservative for liner entry pressure but acceptable.
        m_air_vap = self.res.get('split_primary_vap', 0.0)   # set by zonal_analysis
        A_outer_feed_actual = math.pi * (r_casing_id ** 2 - r_outer_liner_od ** 2)
        A_inner_feed_actual = math.pi * (r_inner_liner_id ** 2 - r_tunnel_od ** 2)
        m_outer_entry = m_air * f_outer_feed              # total entering outer annulus
        m_inner_entry = m_air * (1.0 - f_outer_feed)     # total entering inner annulus
        v_outer_actual = m_outer_entry / (rho2 * A_outer_feed_actual) if A_outer_feed_actual > 0 else 0
        v_inner_actual = m_inner_entry / (rho2 * A_inner_feed_actual) if A_inner_feed_actual > 0 else 0
        # Post-extraction outer annulus velocity (after vaporizer scoops have taken m_air_vap)
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
        self.res['v_outer_post_vap']      = v_outer_post_vap   # after vaporizer extraction
        self.res['outer_liner_circumf_m'] = math.pi * outer_liner_id
        self.res['inner_liner_circumf_m'] = math.pi * inner_liner_od

        # Zone lengths
        self.res['L_primary_mm']   = length * self.DESIGN_PARAMS['frac_primary_zone'] * 1000
        self.res['L_secondary_mm'] = length * self.DESIGN_PARAMS['frac_secondary_zone'] * 1000
        self.res['L_dilution_mm']  = length * self.DESIGN_PARAMS['frac_dilution_zone'] * 1000

    def vaporizer_tubes(self):
        """
        Sizes the reverse-flow vaporizer tubes.

        In the KJ66 architecture:
        - Tubes sit at the REAR dome of the combustion annulus
        - Dome scoops draw both liquid fuel AND a fraction of primary air
          from the outer feed annulus (air-assist vaporization)
        - The combined air+fuel mixture heats in the recirculation zone
        - Exit is bent 180° so the mixture injects FORWARD into the primary zone
        - Counter-injection creates the primary recirculation zone (flame holder)

        Tube count: computed from mean circumference at target pitch, then snapped
        to the nearest value in DESIGN_PARAMS['allowed_vaporizer_counts'] to ensure
        symmetric angular spacing for ignition propagation and thermal balance.

        Exit velocity: uses TOTAL mass flow per tube (fuel + vaporizer air) at the
        corrected exit gas temperature. This is the physically correct calculation —
        prior versions modeled fuel-only, which underestimated exit velocity.
        """
        mean_circumf_m = math.pi * self.res['D_mean_comb_mm'] / 1000
        pitch_m        = self.DESIGN_PARAMS['vap_pitch_mm'] / 1000
        n_raw          = mean_circumf_m / pitch_m

        # Snap to nearest allowed count
        allowed = self.DESIGN_PARAMS['allowed_vaporizer_counts']
        n_tubes = min(allowed, key=lambda n: abs(n - n_raw))

        # --- Vaporizer air co-flow (from zonal_analysis) ---
        m_air_vap       = self.res['split_primary_vap']    # total air via all tubes
        m_fuel_total    = self.res['mdot_fuel']
        m_total_vap     = m_fuel_total + m_air_vap         # combined mass through all tubes

        m_mix_per_tube  = m_total_vap / n_tubes            # total mass per tube
        m_fuel_per_tube = m_fuel_total / n_tubes           # fuel only per tube

        # --- Tube bore sizing on LIQUID FUEL velocity only ---
        # The tube bore carries liquid fuel. Air enters via the dome scoop into the
        # annular gap around the tube body, travelling alongside the vaporizing fuel.
        # Bore diameter is therefore sized on fuel mass flow only.
        rho_liq      = self.FUEL['RHO_LIQ']
        target_v_liq = self.DESIGN_PARAMS['target_tube_liq_vel']
        A_bore       = m_fuel_per_tube / (rho_liq * target_v_liq)
        id_tube      = math.sqrt(4 * A_bore / math.pi)
        id_tube      = max(id_tube, 0.004)   # 4mm minimum (drill / tube availability)
        od_tube      = id_tube + 0.0012      # ~1.2mm wall

        # --- Exit velocity — combined air+fuel through the crimped exit nozzle ---
        # At the 180° bend the air and fuel vapor have merged. The exit nozzle
        # (crimp) must pass the full combined mass. Rather than assuming the exit
        # area equals the bore area, we SIZE the exit nozzle to achieve a target
        # uncrimped velocity, then apply the crimp factor on top.
        #
        # Target uncrimped exit velocity: 40–60 m/s (before 2.5x crimp → 100–150 m/s)
        # This is determined by recirculation zone momentum requirements, not tube bore.
        T_vapor_exit     = 950.0   # K — mixed gas at tube exit
        rho_mix_exit     = self.res['P2_Pa'] / (self.R * T_vapor_exit)
        target_v_exit_uncrimped = 50.0   # m/s — sized to give 125 m/s crimped (mid-range)
        A_exit_nozzle    = m_mix_per_tube / (rho_mix_exit * target_v_exit_uncrimped)
        d_exit_nozzle    = math.sqrt(4 * A_exit_nozzle / math.pi)  # exit nozzle ID

        v_exit_uncrimped = target_v_exit_uncrimped   # by construction
        v_exit_crimped   = v_exit_uncrimped * 2.5

        self.res['vap_n']                = n_tubes
        self.res['vap_n_raw']            = n_raw
        self.res['vap_od_mm']            = od_tube * 1000
        self.res['vap_id_mm']            = id_tube * 1000
        self.res['vap_exit_nozzle_mm']   = d_exit_nozzle * 1000   # exit nozzle ID (crimp point)
        self.res['vap_pitch_actual_mm']  = mean_circumf_m / n_tubes * 1000
        self.res['vap_exit_v_uncrimped'] = v_exit_uncrimped
        self.res['vap_exit_v_crimped']   = v_exit_crimped
        self.res['vap_bend_radius_mm']   = od_tube * 1000 * 1.5
        self.res['vap_m_air_per_tube']   = m_air_vap / n_tubes * 1000
        self.res['vap_m_fuel_per_tube']  = m_fuel_per_tube * 1000
        self.res['vap_AFR_tube']         = (m_air_vap / n_tubes) / m_fuel_per_tube

        # --- Vaporizer Scoop Inlet Sizing (HIGH-SEVERITY — ties CAD geometry to mass model) ---
        # The scoop is the dome opening that draws m_air_vap from the outer feed annulus.
        # Sized using continuity through an orifice with bell-mouth Cd.
        # Output d_scoop_mm is a real fabrication dimension — drill or punch this diameter.
        m_air_vap_total  = self.res.get('split_primary_vap', 0.0)
        m_air_vap_per    = m_air_vap_total / n_tubes if n_tubes > 0 else 0.0
        rho_ann          = self.res['rho2']
        v_scoop_target   = self.DESIGN_PARAMS['vap_scoop_target_vel_m_s']
        Cd_scoop         = self.DESIGN_PARAMS['vap_scoop_cd']
        A_scoop_m2       = m_air_vap_per / (rho_ann * v_scoop_target * Cd_scoop)
        d_scoop_mm       = math.sqrt(4 * A_scoop_m2 / math.pi) * 1000
        # Face velocity at scoop opening = v_scoop_target (by construction — this is what we sized to)
        # The Cd accounts for contraction: effective mass discharge = Cd × A × rho × v_face
        self.res['vap_scoop_dia_mm']      = round(d_scoop_mm, 2)
        self.res['vap_scoop_vel_actual']  = round(v_scoop_target, 1)   # face velocity at scoop inlet

        # --- Fuel pressure drop model ---
        # Two main contributors: (1) straight tube friction, (2) crimp nozzle
        #
        # (1) Straight tube Darcy-Weisbach (liquid fuel, laminar assumed for micro-GT):
        #     Re = rho_liq * v_liq * d_i / mu_fuel
        #     mu_fuel ≈ 2.0e-3 Pa·s  (Jet-A at ~300 K; decreases toward T_boil)
        #     f = 64/Re  (laminar); use Moody for Re > 2300
        #     ΔP_tube = f * (L_tube/d_i) * (rho * v²/2)
        # Effective tube length: straight run ≈ 1.5x chamber length + 2x OD bend
        mu_fuel   = 2.0e-3        # Pa·s — dynamic viscosity, Jet-A at ~300K (conservative)
        Re_tube   = rho_liq * target_v_liq * id_tube / mu_fuel
        if Re_tube < 2300:
            f_darcy = 64.0 / Re_tube   # laminar — common in these tiny tubes
        else:
            # Colebrook approximation (smooth tube, turbulent)
            f_darcy = 0.316 * Re_tube ** (-0.25)   # Blasius

        L_tube_eff = (self.res.get('chamber_length_mm', 120) / 1000) * 1.5 + 2 * od_tube
        dP_tube_Pa = f_darcy * (L_tube_eff / id_tube) * (rho_liq * target_v_liq ** 2 / 2)

        # (2) Crimp nozzle pressure drop:
        #     The crimp creates a sharp-edged contraction. Model as orifice.
        #     Crimp area ratio ≈ 1/2.5 (to achieve 2.5x velocity) → Cd_crimp ≈ 0.61
        #     ΔP_crimp = rho_vapor * v_exit_crimped² / (2 * Cd²)
        Cd_crimp   = 0.61
        dP_crimp   = rho_mix_exit * v_exit_crimped ** 2 / (2 * Cd_crimp ** 2)

        dP_total_Pa  = dP_tube_Pa + dP_crimp
        dP_total_bar = dP_total_Pa / 1e5

        self.res['fuel_Re_tube']      = Re_tube
        self.res['fuel_dP_tube_Pa']   = dP_tube_Pa
        self.res['fuel_dP_crimp_Pa']  = dP_crimp
        self.res['fuel_dP_total_bar'] = dP_total_bar
        self.res['fuel_flow_regime']  = "Laminar" if Re_tube < 2300 else "Turbulent"

    def hole_sizing(self):
        """
        Sizes liner air admission holes for Primary, Secondary, and Dilution zones
        on BOTH the outer liner and inner liner independently.

        METHOD:
        1. Total hole area for each zone from incompressible orifice equation:
               A_total = m_dot / (Cd * sqrt(2 * rho * dP))

        2. Split total area proportionally to liner circumference:
               f_outer = C_outer / (C_outer + C_inner)  (Lefebvre method)
               f_inner = 1 - f_outer

        3. Number of holes on each liner = n_vaporizers × zone_multiplier
               Primary: 2 holes per vaporizer per liner (interspersed between tubes)
               Secondary: 1 hole per vaporizer per liner
               Dilution: 1 hole per vaporizer per liner

        4. Hole diameter from area / count

        5. JET PENETRATION CHECK (Sturgess correlation):
               y_max / H_comb = 1.15 * sqrt(J) * (d_j / H_comb)
               J = (rho_jet * V_jet²) / (rho_gas * V_gas²)
               Target: y_max ≈ 0.5 * H_comb  (reach combustion annulus midplane)

        The penetration check is diagnostic — it tells the engineer if the jets
        will mix thoroughly or just hug the liner wall.
        """
        dP   = self.res['P2_Pa'] * self.DESIGN_PARAMS['target_pressure_drop']
        Cd   = self.DESIGN_PARAMS['discharge_coeff_hole']
        rho2 = self.res['rho2']

        # Flow factor (identical for all holes given same dP, rho, Cd)
        flow_factor = Cd * math.sqrt(2 * rho2 * dP)

        # --- Total hole areas per zone ---
        # Primary holes only carry the liner portion of primary air —
        # the vaporizer fraction is already delivered via the tube scoops.
        A_pri = self.res['split_primary_liner'] / flow_factor
        A_sec = self.res['split_secondary']     / flow_factor
        A_dil = self.res['split_dilution']      / flow_factor

        # --- Circumference-proportional split ---
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

        # --- Hole counts ---
        # Primary: 2 holes per vaporizer pitch per liner
        n_pri_out = n_vap * 2
        n_pri_in  = n_vap * 2

        # Secondary: 1 hole per vaporizer pitch per liner
        n_sec_out = n_vap
        n_sec_in  = n_vap

        # Dilution: start at 2 per vaporizer, scale up if holes would be too large
        # Target max hole dia ≈ 16mm (practical for stainless liner fabrication)
        max_hole_mm = 16.0
        n_dil_out_min = n_vap * 2
        d_dil_out_test = hole_dia(A_dil, f_outer, n_dil_out_min) * 1000
        if d_dil_out_test > max_hole_mm:
            n_dil_out_min = math.ceil(
                (A_dil * f_outer) / (math.pi * (max_hole_mm / 2000) ** 2)
            )
            n_dil_out_min = math.ceil(n_dil_out_min / n_vap) * n_vap

        n_dil_out = n_dil_out_min
        n_dil_in  = n_dil_out

        # --- Hole diameters ---
        d_pri_out = hole_dia(A_pri, f_outer, n_pri_out)
        d_pri_in  = hole_dia(A_pri, f_inner, n_pri_in)
        d_sec_out = hole_dia(A_sec, f_outer, n_sec_out)
        d_sec_in  = hole_dia(A_sec, f_inner, n_sec_in)
        d_dil_out = hole_dia(A_dil, f_outer, n_dil_out)
        d_dil_in  = hole_dia(A_dil, f_inner, n_dil_in)

        # --- Mach number check on hole flow ---
        # Incompressible orifice equation valid for Ma_hole < 0.3.
        # At Ma ≥ 0.3, compressibility reduces discharge area by ~5–10%,
        # meaning true required area is slightly larger than computed here.
        # This is a diagnostic flag; the incompressible result is still used
        # for sizing (conservative for hole area, slightly under-sizes holes).
        T2       = self.res['T2_K']
        a_sound  = math.sqrt(self.GAMMA * self.R * T2)   # speed of sound at inlet conditions

        def hole_mach(m_zone, n_holes, d_hole):
            A_hole = math.pi * (d_hole / 2) ** 2
            if A_hole <= 0 or n_holes <= 0:
                return 0.0
            v_hole = m_zone / (rho2 * A_hole * n_holes)
            return v_hole / a_sound

        # Compute per-hole Mach for each zone (outer liner, as representative)
        # These are stored as diagnostics; no geometry is changed.
        Ma_pri = hole_mach(self.res['split_primary_liner'], n_pri_out, d_pri_out)
        Ma_sec = hole_mach(self.res['split_secondary'], n_sec_out, d_sec_out)
        Ma_dil = hole_mach(self.res['split_dilution'],  n_dil_out, d_dil_out)

        self.res['Ma_hole_pri'] = Ma_pri
        self.res['Ma_hole_sec'] = Ma_sec
        self.res['Ma_hole_dil'] = Ma_dil
        self.res['hole_compressible_warning'] = any(m > 0.3 for m in [Ma_pri, Ma_sec, Ma_dil])

        # --- Explicit Film-Cooling Pattern (HIGH-SEVERITY — mandatory for liner survival) ---
        # 12% reserved air must be translated into actual hole rows with real diameters.
        # Without this, "12% film cooling" exists only as a budget number — the liner
        # has no hardware to protect it from hot-streak impingement.
        #
        # Method: total film area from orifice equation (same dP, same Cd, same rho2),
        # distributed equally across all rows, split by liner circumference fraction.
        # Each row has a fixed hole diameter (1.2 mm standard for micro-jet effusion).
        m_film       = self.res['m_film_cooling']
        A_film_total = m_film / flow_factor   # total effusion area needed

        n_rows_pri = self.DESIGN_PARAMS['film_rows_primary']
        n_rows_sec = self.DESIGN_PARAMS['film_rows_secondary']
        n_rows_dil = self.DESIGN_PARAMS['film_rows_dilution']
        n_rows     = n_rows_pri + n_rows_sec + n_rows_dil

        A_per_row   = A_film_total / max(n_rows, 1)
        d_film_m    = self.DESIGN_PARAMS['film_hole_dia_mm'] / 1000.0
        a_film_hole = math.pi * (d_film_m / 2) ** 2

        # Holes per row on each liner, proportional to circumference (same split as jets)
        n_film_out_per_row = math.ceil(A_per_row * f_outer / a_film_hole)
        n_film_in_per_row  = math.ceil(A_per_row * f_inner / a_film_hole)

        # Axial positions of each cooling row (from combustion zone start = dome)
        pitch_m = self.DESIGN_PARAMS['film_row_pitch_mm'] / 1000.0
        row_positions = [pitch_m * (i + 0.5) for i in range(n_rows)]  # offset from zone boundary

        self.res['film_n_rows']              = n_rows
        self.res['film_n_rows_pri']          = n_rows_pri
        self.res['film_n_rows_sec']          = n_rows_sec
        self.res['film_n_rows_dil']          = n_rows_dil
        self.res['film_hole_dia_mm_actual']  = self.DESIGN_PARAMS['film_hole_dia_mm']
        self.res['film_holes_per_row_outer'] = n_film_out_per_row
        self.res['film_holes_per_row_inner'] = n_film_in_per_row
        self.res['film_total_area_mm2']      = round(A_film_total * 1e6, 1)
        self.res['film_total_holes']         = n_rows * (n_film_out_per_row + n_film_in_per_row)
        # Approximate: after secondary zone, T ≈ 1200 K, P ≈ P2 (approx)
        T_gas_dil  = 1200.0
        rho_gas    = self.res['P2_Pa'] / (self.R * T_gas_dil)

        # Combustion annulus bulk velocity (continuity)
        m_at_dil   = self.res['split_primary'] + self.res['split_secondary'] + self.res['mdot_fuel']
        A_comb     = self.res['combustion_annulus_A']
        cp_dil     = self._get_cp(T_gas_dil)
        # Use rho at dil zone conditions
        V_gas      = m_at_dil / (rho_gas * A_comb)

        # Jet velocity at outer liner dilution holes
        A_hole_out = math.pi * (d_dil_out / 2) ** 2
        V_jet_out  = (self.res['split_dilution'] * f_outer / n_dil_out) / (rho2 * A_hole_out)

        # Jet velocity at inner liner dilution holes
        A_hole_in  = math.pi * (d_dil_in / 2) ** 2
        V_jet_in   = (self.res['split_dilution'] * f_inner / n_dil_in) / (rho2 * A_hole_in)

        # Momentum flux ratios
        J_out = (rho2 * V_jet_out ** 2) / (rho_gas * V_gas ** 2) if V_gas > 0 else 0
        J_in  = (rho2 * V_jet_in  ** 2) / (rho_gas * V_gas ** 2) if V_gas > 0 else 0

        # Penetration depth (Sturgess correlation)
        H_comb = self.res['combustion_gap_mm'] / 1000   # radial height of combustion annulus (m)
        # y_max/H = 1.15 * sqrt(J) * (d_j/H)
        pen_out_norm = 1.15 * math.sqrt(J_out) * (d_dil_out / H_comb) if H_comb > 0 else 0
        pen_in_norm  = 1.15 * math.sqrt(J_in)  * (d_dil_in  / H_comb) if H_comb > 0 else 0
        pen_out_mm   = pen_out_norm * H_comb * 1000
        pen_in_mm    = pen_in_norm  * H_comb * 1000

        # --- Store results ---
        self.res['hole_split_f_outer']    = f_outer
        self.res['hole_split_f_inner']    = f_inner

        # Outer liner holes
        self.res['pri_out_qty']  = n_pri_out
        self.res['pri_out_mm']   = d_pri_out * 1000
        self.res['sec_out_qty']  = n_sec_out
        self.res['sec_out_mm']   = d_sec_out * 1000
        self.res['dil_out_qty']  = n_dil_out
        self.res['dil_out_mm']   = d_dil_out * 1000

        # Inner liner holes
        self.res['pri_in_qty']   = n_pri_in
        self.res['pri_in_mm']    = d_pri_in  * 1000
        self.res['sec_in_qty']   = n_sec_in
        self.res['sec_in_mm']    = d_sec_in  * 1000
        self.res['dil_in_qty']   = n_dil_in
        self.res['dil_in_mm']    = d_dil_in  * 1000

        # Jet penetration diagnostics
        self.res['dil_J_outer']         = J_out
        self.res['dil_J_inner']         = J_in
        self.res['dil_pen_outer_mm']    = pen_out_mm
        self.res['dil_pen_inner_mm']    = pen_in_mm
        self.res['H_comb_mm']           = H_comb * 1000
        self.res['dil_pen_norm_outer']  = pen_out_norm
        self.res['dil_pen_norm_inner']  = pen_in_norm

    def liner_structural(self):
        """
        Two structural checks for the thin combustor liner walls:

        A) HOOP STRESS (pressure vessel, thin-wall assumption):
               σ_hoop = ΔP × r / t
           where r = liner mean radius, t = wall thickness, ΔP = pressure drop across liner.
           For 304 SS at 900°C: yield strength ≈ 90–130 MPa.
           For Inconel 625 at 900°C: yield strength ≈ 350–480 MPa.
           Safety factor of 2x applied against yield.

        B) MINIMUM LIGAMENT SPACING:
               e_min = (circumference / n_holes) - d_hole
           Minimum recommended: 2× wall thickness or 1.5mm, whichever is larger.
           Below this, stress concentration between adjacent holes risks tearing
           during thermal cycling (a common failure mode in micro-jet liners).

        NOTE: These are first-order thin-wall checks. They do not account for:
           - Thermal gradient stresses (significant in thin stainless)
           - Creep at sustained high temperature
           - Weld heat affected zones (if liner is welded / rolled)
        Treat results as screening — if flagged, perform detailed FEA.
        """
        wall_m = self.inputs['wall_thickness_mm'] / 1000.0
        dP     = self.res['P2_Pa'] * self.DESIGN_PARAMS['target_pressure_drop']

        # --- Hoop stress: outer liner (larger radius = higher stress) ---
        r_outer_mean = (self.res['outer_liner_od_mm'] + self.res['outer_liner_id_mm']) / 2 / 1000 / 2
        sigma_outer  = dP * r_outer_mean / wall_m / 1e6   # MPa

        # --- Hoop stress: inner liner ---
        r_inner_mean = (self.res['inner_liner_od_mm'] + self.res['inner_liner_id_mm']) / 2 / 1000 / 2
        sigma_inner  = dP * r_inner_mean / wall_m / 1e6   # MPa

        # Material limits (conservative, at operating temperature ~900°C)
        sigma_allow_ss304   =  65.0  # MPa — 304 SS hot yield / 2 (SF=2)
        sigma_allow_in625   = 175.0  # MPa — Inconel 625 hot yield / 2 (SF=2)

        # --- Minimum ligament spacing ---
        # Check dilution holes (largest diameter, most critical)
        # Outer liner
        circ_outer_mm = self.res['outer_liner_circumf_m'] * 1000
        pitch_out_mm  = circ_outer_mm / self.res['dil_out_qty']
        lig_out_mm    = pitch_out_mm - self.res['dil_out_mm']

        # Inner liner
        circ_inner_mm = self.res['inner_liner_circumf_m'] * 1000
        pitch_in_mm   = circ_inner_mm / self.res['dil_in_qty']
        lig_in_mm     = pitch_in_mm - self.res['dil_in_mm']

        lig_min_req   = max(self.inputs['wall_thickness_mm'] * 2.0, 1.5)   # mm

        self.res['sigma_hoop_outer_MPa']  = sigma_outer
        self.res['sigma_hoop_inner_MPa']  = sigma_inner
        self.res['sigma_allow_ss304_MPa'] = sigma_allow_ss304
        self.res['sigma_allow_in625_MPa'] = sigma_allow_in625
        self.res['ligament_outer_mm']     = lig_out_mm
        self.res['ligament_inner_mm']     = lig_in_mm
        self.res['ligament_min_req_mm']   = lig_min_req
        self.res['liner_ss304_ok']        = (sigma_outer < sigma_allow_ss304 and
                                              sigma_inner < sigma_allow_ss304)
        self.res['liner_in625_ok']        = (sigma_outer < sigma_allow_in625 and
                                              sigma_inner < sigma_allow_in625)
        self.res['ligament_outer_ok']     = lig_out_mm >= lig_min_req
        self.res['ligament_inner_ok']     = lig_in_mm  >= lig_min_req


    def temperature_traverse_quality(self):
        """
        ⚠ MIXING INDEX — NOT A TRUE TTQ PREDICTION ⚠

        True Temperature Traverse Quality (TTQ) requires:
          - CFD or rig-measured temperature profiles at the turbine inlet plane
          - Swirl number, circumferential hole spacing, and dome geometry inputs
          - Dilution jet-swirl interaction modeling

        This function computes a first-order MIXING INDEX (MI) based solely on
        dilution jet penetration depth relative to the combustion annulus midplane.
        It is a qualitative screening metric only.

            MI target < 0.25 (analogous to TTQ design target for annular combustors)

        Perfect mixing (MI → 0.05): jets from both walls reach exactly midplane.
        Poor mixing   (MI → 0.40): jets barely penetrate, hot core survives to TIT.
        """
        pen_norm_out = self.res['dil_pen_norm_outer']
        pen_norm_in  = self.res['dil_pen_norm_inner']

        eff_out = 1.0 - abs(pen_norm_out - 0.5) / 0.5
        eff_in  = 1.0 - abs(pen_norm_in  - 0.5) / 0.5
        eff_out = max(0.0, min(1.0, eff_out))
        eff_in  = max(0.0, min(1.0, eff_in))
        mix_eff = (eff_out + eff_in) / 2.0

        MI = 0.40 - mix_eff * 0.35
        self.res['mixing_index']  = MI
        self.res['mix_eff_outer'] = eff_out
        self.res['mix_eff_inner'] = eff_in
    def get_cad_geometry(self):
        """
        Returns ONE clean dictionary with ALL geometry for CAD import.
        Use like:
            cad = model.get_cad_geometry()
            json.dump(cad, open('combustor_cad.json', 'w'), indent=2)
        All linear dimensions in mm. Counts and fractions are unitless.
        """
        return {
            "metadata": {
                "units": "mm",
                "description": "Reverse-flow annular micro-jet combustor — KJ66 style",
                "generated_by": "V7_CombustionChamberDesign.py"
            },
            "overall": {
                "chamber_length": self.res['chamber_length_mm'],
                "mean_comb_dia":  self.res['D_mean_comb_mm']
            },
            "casing": {
                "od": self.res['casing_od_mm'],
                "id": self.res['casing_id_mm'],
                "wall_thickness": self.inputs['wall_thickness_mm']
            },
            "shaft_tunnel": {
                "od": self.res['shaft_tunnel_od_mm'],
                "id": self.res['shaft_tunnel_id_mm']
            },
            "outer_liner": {
                "od": self.res['outer_liner_od_mm'],
                "id": self.res['outer_liner_id_mm'],
                "length": self.res['chamber_length_mm'],
                "annulus_gap": self.res['outer_annulus_gap_mm']
            },
            "inner_liner": {
                "od": self.res['inner_liner_od_mm'],
                "id": self.res['inner_liner_id_mm'],
                "length": self.res['chamber_length_mm'],
                "annulus_gap": self.res['inner_annulus_gap_mm']
            },
            "combustion_annulus": {
                "radial_height": self.res['combustion_gap_mm'],
                "mean_diameter": self.res['D_mean_comb_mm'],
                "cross_section_area_mm2": self.res['combustion_annulus_A'] * 1e6
            },
            "vaporizers": {
                "count": self.res['vap_n'],
                "raw_count": self.res['vap_n_raw'],
                "pitch_actual": self.res['vap_pitch_actual_mm'],
                "tube_od": self.res['vap_od_mm'],
                "tube_bore_id": self.res['vap_id_mm'],
                "exit_nozzle_id": self.res['vap_exit_nozzle_mm'],
                "scoop_inlet_dia": self.res['vap_scoop_dia_mm'],
                "bend_radius": self.res['vap_bend_radius_mm'],
                "per_tube": {
                    "fuel_g_s": self.res['vap_m_fuel_per_tube'],
                    "air_g_s": self.res['vap_m_air_per_tube']
                }
            },
            "main_holes": {
                "outer_fraction": self.res['hole_split_f_outer'],
                "inner_fraction": self.res['hole_split_f_inner'],
                "primary": {
                    "outer_qty": self.res['pri_out_qty'],
                    "outer_dia": self.res['pri_out_mm'],
                    "inner_qty": self.res['pri_in_qty'],
                    "inner_dia": self.res['pri_in_mm']
                },
                "secondary": {
                    "outer_qty": self.res['sec_out_qty'],
                    "outer_dia": self.res['sec_out_mm'],
                    "inner_qty": self.res['sec_in_qty'],
                    "inner_dia": self.res['sec_in_mm']
                },
                "dilution": {
                    "outer_qty": self.res['dil_out_qty'],
                    "outer_dia": self.res['dil_out_mm'],
                    "inner_qty": self.res['dil_in_qty'],
                    "inner_dia": self.res['dil_in_mm']
                }
            },
            "film_cooling": {
                "total_rows": self.res['film_n_rows'],
                "rows_primary": self.res['film_n_rows_pri'],
                "rows_secondary": self.res['film_n_rows_sec'],
                "rows_dilution": self.res['film_n_rows_dil'],
                "row_pitch": self.DESIGN_PARAMS['film_row_pitch_mm'],
                "hole_dia": self.res['film_hole_dia_mm_actual'],
                "outer_holes_per_row": self.res['film_holes_per_row_outer'],
                "inner_holes_per_row": self.res['film_holes_per_row_inner'],
                "total_holes": self.res['film_total_holes'],
                "total_area_mm2": self.res['film_total_area_mm2']
            },
            "annulus_velocities": {
                "outer_entry": self.res['v_outer_annulus'],
                "outer_post_vap": self.res['v_outer_post_vap'],
                "inner": self.res['v_inner_annulus']
            },
            "zone_lengths": {
                "primary": self.res['L_primary_mm'],
                "secondary": self.res['L_secondary_mm'],
                "dilution": self.res['L_dilution_mm']
            }
        }
    def run(self):
        self.thermodynamics()
        self.mass_flow_and_fuel()
        self.zonal_analysis()
        self.mechanical_geometry()
        self.vaporizer_tubes()
        self.hole_sizing()
        self.liner_structural()
        self.temperature_traverse_quality()
        self.res['cad_geometry'] = self.get_cad_geometry()
        return self.res
        return self.res

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
    print("   KJ66-Style Architecture  |  Rev 10 — Physics Corrections")
    print("=" * W)

    # -- ASCII cross-section --
    oc_mm  = res['casing_od_mm']
    ol_mm  = res['outer_liner_od_mm']
    il_mm  = res['inner_liner_od_mm']
    st_mm  = res['shaft_tunnel_od_mm']
    scale  = oc_mm / 30.0  # chars per mm, ~30 chars wide for outer casing

    def bar(od, id_=None):
        """Simple 1-D cross section bar, normalized to outer casing."""
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

    # [1] Thermodynamics
    header("[1] THERMODYNAMICS")
    row("Inlet Total Pressure P2", f"{res['P2_Pa']/1000:.1f}", "kPa")
    row("Inlet Total Temperature T2", f"{res['T2_K']:.0f}", "K")
    row("Target Turbine Inlet Temp", f"{inputs['target_tit_k']:.0f}", "K")
    clamp_note = "CLAMPED — real value higher" if res['T_primary_zone_clamped'] else "rich φ=1.6"
    row("Est. Primary Zone Temperature", f"{res['T_primary_zone_est']:.0f}", "K", clamp_note)
    if res['T_primary_zone_clamped']:
        print(f"  (Raw calculated: {res['T_primary_zone_raw']:.0f} K — clamped to 2050 K physical limit)")
    row("Avg Cp Used", f"{res['cp_used']:.1f}", "J/kg·K")
    row("Inlet Air Density ρ2", f"{res['rho2']:.3f}", "kg/m³")

    # [2] Flow Budget
    header("[2] FLOW BUDGET")
    m_total = res['mdot_air']
    row("Total Air Mass Flow", f"{m_total:.4f}", "kg/s")
    row("Total Fuel Mass Flow", f"{res['mdot_fuel']*1000:.2f}", "g/s")
    row("Overall AFR", f"{res['overall_AFR']:.1f}", "kg_air/kg_fuel")
    row("Overall Equivalence Ratio φ", f"{res['overall_phi']:.3f}", "")
    print()
    print("  Air Split:")
    row("  Film Cooling (Effusion)",
        f"{res['m_film_cooling']/m_total:.1%}", "",
        f"{res['m_film_cooling']*1000:.1f} g/s  ({res['film_n_rows']} rows)")
    row("  Primary — Vaporizer Scoops",
        f"{res.get('split_primary_vap',0)/m_total:.1%}", "",
        f"{res.get('split_primary_vap',0)*1000:.1f} g/s")
    row("  Primary — Liner Holes",
        f"{res['split_primary_liner']/m_total:.1%}", "",
        f"{res['split_primary_liner']*1000:.1f} g/s")
    row("  Secondary Zone",
        f"{res['split_secondary']/m_total:.1%}", "",
        f"{res['split_secondary']*1000:.1f} g/s  φ→0.60")
    row("  Dilution Zone",
        f"{res['split_dilution']/m_total:.1%}", "",
        f"{res['split_dilution']*1000:.1f} g/s")

    flag(res['split_dilution'] <= 0,
         "ZERO dilution air — TIT target may exceed what this air mass can achieve "
         "at this equivalence ratio. Check overall φ < 0.6 after secondary zone.")

    # [3] Four-Wall Annular Geometry
    header("[3] FOUR-WALL ANNULAR GEOMETRY")
    print("  OUTER BOUNDARY:")
    row("  Outer Casing OD", f"{res['casing_od_mm']:.2f}", "mm")
    row("  Outer Casing ID", f"{res['casing_id_mm']:.2f}", "mm")
    row("  Outer Feed Annulus Gap", f"{res['outer_annulus_gap_mm']:.2f}", "mm")
    row("  Outer Annulus Velocity", f"{res['v_outer_annulus']:.1f}", "m/s")
    print()
    print("  OUTER LINER:")
    row("  Outer Liner OD", f"{res['outer_liner_od_mm']:.2f}", "mm")
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
    row("  Inner Liner OD", f"{res['inner_liner_od_mm']:.2f}", "mm")
    row("  Inner Liner ID", f"{res['inner_liner_id_mm']:.2f}", "mm")
    print()
    print("  INNER BOUNDARY:")
    row("  Shaft Tunnel OD", f"{res['shaft_tunnel_od_mm']:.2f}", "mm")
    row("  Shaft Tunnel ID", f"{res['shaft_tunnel_id_mm']:.2f}", "mm")
    row("  Inner Feed Annulus Gap", f"{res['inner_annulus_gap_mm']:.2f}", "mm")
    row("  Inner Annulus Velocity", f"{res['v_inner_annulus']:.1f}", "m/s")

    flag(res['v_outer_annulus'] > 50,
         f"Outer annulus velocity {res['v_outer_annulus']:.1f} m/s is high — "
         "risk of pressure distortion. Consider increasing casing OD.")
    flag(res['v_inner_annulus'] > 50,
         f"Inner annulus velocity {res['v_inner_annulus']:.1f} m/s is high — "
         "consider increasing shaft tunnel clearance or reducing inner feed fraction.")
    flag(res['combustion_gap_mm'] < 8,
         f"Combustion annulus height {res['combustion_gap_mm']:.1f} mm is very narrow. "
         "Flame stability and jet penetration may be poor.")

    f_split = res.get('f_outer_feed_used', 0.60)
    print(f"\n  Feed split assumption: {f_split:.0%} outer / {1-f_split:.0%} inner")
    print(f"  (Mass split imposed = area split — see DESIGN_PARAMS['f_outer_feed'])")
    print(f"  Both annuli velocity = {res['v_outer_annulus']:.1f} m/s by construction.")

    # [4] Vaporizer Tubes
    header("[4] VAPORIZER TUBES  (\"Candy Canes\"  — Air-Assist)")
    row("Number of Tubes (snapped)", f"{res['vap_n']}", "",
        f"raw={res['vap_n_raw']:.1f}, from {[6,8,12]}")
    row("Scoop Inlet Diameter", f"{res['vap_scoop_dia_mm']:.2f}", "mm",
        f"@ {res['vap_scoop_vel_actual']:.1f} m/s  (bell-mouth, 15–25° flare)")
    row("Tube OD / Bore ID", f"{res['vap_od_mm']:.2f} / {res['vap_id_mm']:.2f}", "mm")
    row("Exit Nozzle ID (at crimp)", f"{res['vap_exit_nozzle_mm']:.2f}", "mm")
    row("U-bend Radius (1.5×OD)", f"{res['vap_bend_radius_mm']:.1f}", "mm")
    print()
    print("  Per-tube mass flow (air-assist):")
    row("  Fuel per tube", f"{res['vap_m_fuel_per_tube']:.2f}", "g/s")
    row("  Air per tube (scoop)", f"{res['vap_m_air_per_tube']:.2f}", "g/s")
    row("  Tube AFR", f"{res['vap_AFR_tube']:.2f}", "", "φ_tube from primary air fraction")
    print()
    row("Exit Vel (uncrimped, air+fuel)", f"{res['vap_exit_v_uncrimped']:.1f}", "m/s",
        "at T_exit=950K, total mass")
    row("Exit Vel (crimped ×2.5)", f"{res['vap_exit_v_crimped']:.1f}", "m/s",
        "target 50–120 m/s")
    print()
    print("  Outer annulus velocity accounting:")
    row("  Entry velocity (pre-scoop)", f"{res['v_outer_annulus']:.1f}", "m/s")
    row("  Post-scoop velocity", f"{res['v_outer_post_vap']:.1f}", "m/s",
        "after vaporizer extraction")
    print()
    print("  Fuel Pressure Drop (per tube):")
    row("  Flow Regime", f"  {res['fuel_flow_regime']}", "",
        f"Re = {res['fuel_Re_tube']:.0f}")
    row("  Tube friction ΔP", f"{res['fuel_dP_tube_Pa']/1e5:.4f}", "bar")
    row("  Crimp nozzle ΔP", f"{res['fuel_dP_crimp_Pa']/1e5:.4f}", "bar")
    row("  TOTAL fuel ΔP", f"{res['fuel_dP_total_bar']:.3f}", "bar",
        "typical KJ66 range: 0.4–0.8 bar")

    flag(res['vap_exit_v_crimped'] < 50,
         "Crimped exit velocity < 50 m/s — flashback risk in the vaporizer bend.")
    flag(res['vap_exit_v_crimped'] > 150,
         "Crimped exit velocity > 150 m/s — may approach choked flow; reduce crimp ratio.")
    flag(res['vap_id_mm'] <= 4.1,
         "Tube bore at minimum (4 mm). Consider lower liquid fuel velocity target.")
    flag(res['fuel_dP_total_bar'] < 0.3,
         f"Fuel ΔP = {res['fuel_dP_total_bar']:.3f} bar may be too low for stable "
         "atomization — fuel pump head may be marginal.")
    flag(res['fuel_dP_total_bar'] > 1.2,
         f"Fuel ΔP = {res['fuel_dP_total_bar']:.3f} bar is high — verify pump selection "
         "and check tube bore isn't undersized.")

    print()
    print("  NOTE: Vaporizer tubes feed from the OUTER feed annulus at the rear dome.")
    print("  The 180° bend injects fuel FORWARD (counter-flow) to create the")
    print("  primary recirculation zone (aerodynamic flame holder).")

    # [5] Hole Sizing — DUAL LINER
    header("[5] LINER HOLE SIZING  (Dual-Wall Annular)")
    print(f"  Hole area split — Outer: {res['hole_split_f_outer']:.1%}  Inner: {res['hole_split_f_inner']:.1%}")
    print(f"  (Proportional to liner circumference — Lefebvre method)")
    print()

    def hole_row(zone, qty_out, d_out, qty_in, d_in):
        print(f"  {zone:12s}  OUTER: {qty_out:2d} × {d_out:5.2f} mm     INNER: {qty_in:2d} × {d_in:5.2f} mm")

    hole_row("PRIMARY:",
             res['pri_out_qty'], res['pri_out_mm'],
             res['pri_in_qty'],  res['pri_in_mm'])
    hole_row("SECONDARY:",
             res['sec_out_qty'], res['sec_out_mm'],
             res['sec_in_qty'],  res['sec_in_mm'])
    hole_row("DILUTION:",
             res['dil_out_qty'], res['dil_out_mm'],
             res['dil_in_qty'],  res['dil_in_mm'])

    total_holes = (res['pri_out_qty'] + res['pri_in_qty'] +
                   res['sec_out_qty'] + res['sec_in_qty'] +
                   res['dil_out_qty'] + res['dil_in_qty'])
    print(f"\n  Total primary/secondary/dilution holes (both liners): {total_holes}")

    print()
    print("  Hole Mach Numbers (outer liner, incompressible eq. valid for Ma < 0.3):")
    print(f"  Pri: Ma={res['Ma_hole_pri']:.3f}  Sec: Ma={res['Ma_hole_sec']:.3f}  "
          f"Dil: Ma={res['Ma_hole_dil']:.3f}")
    flag(res['hole_compressible_warning'],
         "One or more hole zones have Ma > 0.3 — incompressible area equation "
         "underestimates required hole area by ~5–10%. Consider upsizing by that margin "
         "or use compressible orifice equation for final design.")

    # [5b] Explicit Film Cooling Pattern
    header("[5b] LINER FILM COOLING PATTERN  (Effusion Holes)")
    print(f"  {res['film_n_rows']} total rows: "
          f"{res['film_n_rows_pri']} primary / "
          f"{res['film_n_rows_sec']} secondary / "
          f"{res['film_n_rows_dil']} dilution")
    print(f"  Row pitch: {res['film_hole_dia_mm_actual']:.2f} mm holes "
          f"at {inputs.get('film_row_pitch_mm', 38.0):.0f} mm axial spacing")
    print()
    row("Outer liner — holes per row", f"{res['film_holes_per_row_outer']}", "holes",
        f"× {res['film_n_rows']} rows = {res['film_holes_per_row_outer']*res['film_n_rows']} total")
    row("Inner liner — holes per row", f"{res['film_holes_per_row_inner']}", "holes",
        f"× {res['film_n_rows']} rows = {res['film_holes_per_row_inner']*res['film_n_rows']} total")
    row("Total effusion holes (both liners)", f"{res['film_total_holes']}", "")
    row("Total film area", f"{res['film_total_area_mm2']:.1f}", "mm²")
    print()
    print("  Placement rule: first row immediately downstream of each primary jet hole.")
    print("  Stagger adjacent rows by 30° circumferentially to avoid hot streaks.")
    flag(res['film_total_area_mm2'] < 100,
         "Film area very small — verify film_cooling_fraction is not under-reserved.")
    flag(res['film_holes_per_row_outer'] < 8,
         "Fewer than 8 holes per row on outer liner — coverage may be insufficient. "
         "Reduce film_hole_dia_mm or increase film_cooling_fraction.")

    # [6] Jet Penetration Analysis
    header("[6] JET PENETRATION ANALYSIS  (Dilution Zone)")
    print("  Target: jets from each wall reach combustion annulus MIDPLANE")
    print(f"  Combustion annulus height H = {res['H_comb_mm']:.1f} mm")
    print()
    row("  Outer liner momentum flux J", f"{res['dil_J_outer']:.2f}", "")
    row("  Outer jet penetration depth", f"{res['dil_pen_outer_mm']:.1f}", "mm",
        f"{res['dil_pen_norm_outer']:.2f}×H  (ideal: 0.50×H)")
    print()
    row("  Inner liner momentum flux J", f"{res['dil_J_inner']:.2f}", "")
    row("  Inner jet penetration depth", f"{res['dil_pen_inner_mm']:.1f}", "mm",
        f"{res['dil_pen_norm_inner']:.2f}×H  (ideal: 0.50×H)")

    flag(res['dil_pen_norm_outer'] < 0.35,
         f"Outer dil. jets only penetrate {res['dil_pen_norm_outer']:.2f}×H — "
         "increase hole size or reduce n_holes on outer liner.")
    flag(res['dil_pen_norm_outer'] > 0.70,
         f"Outer dil. jets over-penetrate {res['dil_pen_norm_outer']:.2f}×H — "
         "may impinge on inner liner. Reduce hole size or increase n_holes.")
    flag(res['dil_pen_norm_inner'] < 0.35,
         f"Inner dil. jets only penetrate {res['dil_pen_norm_inner']:.2f}×H — "
         "increase hole size or reduce n_holes on inner liner.")

    # [7] Structural Checks
    header("[7] LINER STRUCTURAL CHECKS")
    print("  Hoop stress (thin-wall: σ = ΔP·r/t), SF=2 applied to hot yield strength:")
    row("  Outer liner hoop stress", f"{res['sigma_hoop_outer_MPa']:.1f}", "MPa")
    row("  Inner liner hoop stress", f"{res['sigma_hoop_inner_MPa']:.1f}", "MPa")
    row("  Allowable — 304 SS @ 900°C", f"{res['sigma_allow_ss304_MPa']:.0f}", "MPa")
    row("  Allowable — Inconel 625 @ 900°C", f"{res['sigma_allow_in625_MPa']:.0f}", "MPa")

    ok_ss = "✓ OK" if res['liner_ss304_ok'] else "✗ EXCEEDS"
    ok_in = "✓ OK" if res['liner_in625_ok'] else "✗ EXCEEDS"
    print(f"  304 SS:       {ok_ss}  |  Inconel 625:  {ok_in}")

    flag(not res['liner_in625_ok'],
         "Hoop stress exceeds Inconel 625 allowable — increase wall thickness or "
         "reduce pressure drop target.")
    flag(not res['liner_ss304_ok'] and res['liner_in625_ok'],
         "Hoop stress exceeds 304 SS allowable but within Inconel 625 limits. "
         "Use Inconel 625 (or equivalent high-temp alloy) for liner fabrication.")

    print()
    print("  Minimum ligament between dilution holes (most critical — largest holes):")
    row("  Outer liner ligament", f"{res['ligament_outer_mm']:.1f}", "mm")
    row("  Inner liner ligament", f"{res['ligament_inner_mm']:.1f}", "mm")
    row("  Required minimum", f"{res['ligament_min_req_mm']:.1f}", "mm",
        "= max(2×wall, 1.5mm)")

    flag(not res['ligament_outer_ok'],
         f"Outer liner dilution holes too close — ligament only {res['ligament_outer_mm']:.1f} mm. "
         f"Increase hole count to reduce diameter, or reduce target dP.")
    flag(not res['ligament_inner_ok'],
         f"Inner liner dilution holes too close — ligament only {res['ligament_inner_mm']:.1f} mm. "
         f"Same corrective action as outer liner applies.")
    if res['ligament_outer_ok'] and res['ligament_inner_ok']:
        print("  ✓ Ligament spacing acceptable on both liners.")

    print()
    print("  NOTE: Thermal growth not modeled. Add 0.1–0.3mm radial clearance")
    print("  on liner ODs for differential expansion at operating temperature.")

    # [8] Mixing Index
    header("[8] MIXING INDEX  (first-order, NOT a true TTQ prediction)")
    print("  ⚠ See method docstring — this is a penetration-based screening metric.")
    print("  True TTQ requires CFD or rig traverse data.")
    print()
    row("Mixing Index (MI)", f"{res['mixing_index']:.3f}", "", "target < 0.25")
    row("Outer jet mix effectiveness", f"{res['mix_eff_outer']:.1%}", "")
    row("Inner jet mix effectiveness", f"{res['mix_eff_inner']:.1%}", "")
    flag(res['mixing_index'] > 0.25,
         f"MI = {res['mixing_index']:.3f} > 0.25. Improve dilution jet penetration "
         "to reduce turbine blade hot-spot risk.")

    print()
    print("=" * W)
    print("  DESIGN NOTES:")
    print("  • Liner wall = 0.5–1.5mm stainless steel (304 SS) or Inconel 625")
    print("  • All jet holes: deburred, edge-broken 0.1–0.2mm chamfer")
    print("  • Vaporizer tubes: 316 SS or Inconel 625 for thermal cycling resistance")
    print("  • Vaporizer scoops: 15–25° bell-mouth flare at dome inlet (not sharp-edged)")
    print("  • Film cooling rows: first row immediately after primary holes; stagger 30°")
    print("  • Install mid-span support ring on vaporizer tubes (0.8–1.0mm Inconel wire)")
    print("    — prevents fatigue cracking from combustion pressure oscillations")
    print("  • Inner-liner struts: 3× teardrop cross-section, ≤7% flow blockage,")
    print("    0.20mm cold radial clearance to shaft tunnel for thermal growth")
    print("  • Igniter plug: outer liner, primary zone, midway between vaporizers")
    print("  • Thermal growth: add 0.1–0.3mm radial clearance on liner ODs")
    print("=" * W)
    # [9] CAD Geometry Export
    header("[9] CAD GEOMETRY EXPORT  (copy-paste into CAD)")
    print("  → Use model.get_cad_geometry() in your script")
    print("  → Or copy the JSON below for parametric import")
    print()
    import json
    cad = res['cad_geometry']
    print(json.dumps(cad, indent=2))
    print()


# ==============================================================================
# PRESET CONFIGURATIONS
# ==============================================================================
# Note: shaft_tunnel_od_inch is the OD of the bearing/shaft housing tube
# For KJ66 (110mm engine): shaft tunnel ≈ 30mm OD  →  30/25.4 ≈ 1.18"
# For a 6" (152mm) engine: shaft tunnel scales proportionally ≈ 1.65"

user_inputs_6in = {
    'casing_od_inch':        6.00,   # 152mm
    'shaft_tunnel_od_inch':  1.65,   # ~42mm (scaled from KJ66 ratio)
    'wall_thickness_mm':     1.50,
    'pressure_ratio':        1.50,
    'compressor_efficiency': 0.94,
    'mass_flow_air_kg_s':    0.487,
    'target_tit_k':          900.0,
}

kj66_inputs = {
    'casing_od_inch':        4.33,   # 110mm — KJ66 nominal
    'shaft_tunnel_od_inch':  1.18,   # ~30mm — KJ66 bearing housing OD
    'wall_thickness_mm':     0.50,
    'pressure_ratio':        2.20,
    'compressor_efficiency': 0.74,
    'mass_flow_air_kg_s':    0.23,
    'target_tit_k':         1123.0,
}


# ==============================================================================
# MAIN
# ==============================================================================
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
    print("─" * 65)
    print("  REVERSE-FLOW ANNULAR COMBUSTOR DESIGN TOOL")
    print("  KJ66-style Micro-Jet  |  Rev 10")
    print("─" * 65)
    print()
    print("  Presets available:")
    print("    1 — 6\" custom engine  (PR 1.5,  0.487 kg/s,  TIT 900 K)")
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
        print("  (Enter 'none' for mass_flow to use casing-area autoscale)")
        print()
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