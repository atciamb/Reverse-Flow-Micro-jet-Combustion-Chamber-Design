"""
────────────────────────────────────────────────────────────────────────────
NX expression name               V21 result key / source
────────────────────────────────────────────────────────────────────────────
casing_id                        res['casing_id_mm']
casing_od                        res['casing_od_mm']
chamber_length                   res['chamber_length_mm']
combustion_annulus_radial_height res['combustion_gap_mm']
wall_thickness                   inputs['wall_thickness_mm']
prim_sec_length                  res['L_primary_mm']          ← primary zone only
dil_length                       res['L_dilution_mm']
shaft_tunnel_id                  res['shaft_tunnel_id_mm']
shaft_tunnel_od                  res['shaft_tunnel_od_mm']
outer_liner_id                   res['outer_liner_id_mm']
outer_liner_od                   res['outer_liner_od_cold_mm']  ← COLD BUILD dimension
outer_annulus_gap                res['outer_annulus_gap_mm']
inner_liner_id                   res['inner_liner_id_mm']
inner_liner_od                   res['inner_liner_od_cold_mm']  ← COLD BUILD dimension
inner_liner_annulus_gap          res['inner_annulus_gap_mm']
outer_primHole_dia               res['pri_out_mm']
inner_primHole_dia               res['pri_in_mm']
outer_secHole_dia                res['sec_out_mm']
inner_secHole_dia                res['sec_in_mm']
outer_dilHole_dia                res['dil_out_mm']
inner_dilHole_dia                res['dil_in_mm']
primHole_qty                     res['pri_out_qty']   ← same count on both liners
secHole_qty                      res['sec_out_qty']   ← same count on both liners
dilHole_qty                      res['dil_out_qty']   ← same count on both liners
outer_film_hole_qty              res['film_holes_per_row_outer']  ← per row
inner_film_hole_qty              res['film_holes_per_row_inner']  ← per row
film_hole_area                   res['film_total_area_mm2']
"""

import os
import sys
import importlib.util
from datetime import datetime


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
V21_CANDIDATES = [
    os.path.join(SCRIPT_DIR, 'V21_CombustionChamberDesign.py'),
    os.path.join(SCRIPT_DIR, '..', 'V21_CombustionChamberDesign.py'),
]

v21_path = None
for candidate in V21_CANDIDATES:
    if os.path.isfile(candidate):
        v21_path = os.path.abspath(candidate)
        break

if v21_path is None:
    print("ERROR: Cannot find V21_CombustionChamberDesign.py")
    print("  Looked in:")
    for c in V21_CANDIDATES:
        print(f"    {os.path.abspath(c)}")
    print("  Place V21_CombustionChamberDesign.py in the same folder as this script.")
    sys.exit(1)

# Load the module without running its __main__ block
spec = importlib.util.spec_from_file_location("V21_combustor", v21_path)
v21  = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v21)

MicroJetCombustor = v21.MicroJetCombustor
user_inputs_6in   = v21.user_inputs_6in
kj66_inputs       = v21.kj66_inputs
KEYS              = v21.KEYS


print()
print("─" * 65)
print("  NX JOURNAL GENERATOR  —  V21 Combustor Design Tool")
print("─" * 65)
print()
print("  Presets:")
print("    1 — 6\" custom engine  (PR 1.5,  0.487 kg/s,  TIT 900 K)")
print("    2 — KJ66 reference    (PR 2.2,  0.230 kg/s, TIT 1123 K)")
print("    3 — Manual input")
print()
choice = input("  Select [1/2/3]: ").strip()

if choice == "1":
    chosen_inputs = user_inputs_6in
    preset_label  = "6-inch engine preset"
elif choice == "2":
    chosen_inputs = kj66_inputs
    preset_label  = "KJ66 reference preset"
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
    preset_label = "manual input"

print(f"\n  Inputs: {chosen_inputs}")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Run the combustor model
# ─────────────────────────────────────────────────────────────────────────────

try:
    model = MicroJetCombustor(chosen_inputs)
    res   = model.run()
except ValueError as e:
    print(f"\n  GEOMETRY ERROR: {e}")
    print("  Adjust casing OD, shaft tunnel OD, or wall thickness and retry.")
    sys.exit(1)

print(f"\n  Model run complete.")
print(f"    mdot_fuel : {res['mdot_fuel']*1000:.2f} g/s")
print(f"    phi       : {res['overall_phi']:.3f}")
print(f"    gap       : {res['combustion_gap_mm']:.1f} mm")
print(f"    n_vap     : {res['vap_n']}")

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Build the expression table
#
#     Each entry:  nx_name → (value, unit_tag, comment)
#     unit_tag:
#       'mm'   → workPart.UnitCollection.FindObject("MilliMeter")
#       'mm2'  → workPart.UnitCollection.FindObject("SquareMilliMeter")
#       'int'  → NXOpen.Unit.Null   (dimensionless integer)
# ─────────────────────────────────────────────────────────────────────────────

def fmt_mm(v, dp=3):
    """Round to dp decimal places, return as string."""
    return f"{round(float(v), dp):.{dp}f}"

def fmt_int(v):
    """Integer quantity — no decimal point."""
    return str(int(round(float(v))))

def fmt_mm2(v, dp=1):
    return f"{round(float(v), dp):.{dp}f}"

expressions = [
    # ── section: global casing ───────────────────────────────────────────────
    ("casing_od",
     fmt_mm(res['casing_od_mm']), "mm",
     "Outer casing outside diameter"),

    ("casing_id",
     fmt_mm(res['casing_id_mm']), "mm",
     "Outer casing inside diameter  = casing_od - 2*wall"),

    ("wall_thickness",
     fmt_mm(chosen_inputs['wall_thickness_mm']), "mm",
     "Casing and liner nominal wall thickness"),

    # ── section: shaft tunnel ────────────────────────────────────────────────
    ("shaft_tunnel_od",
     fmt_mm(res['shaft_tunnel_od_mm']), "mm",
     "Shaft tunnel outside diameter"),

    ("shaft_tunnel_id",
     fmt_mm(res['shaft_tunnel_id_mm']), "mm",
     "Shaft tunnel inside diameter  = shaft_tunnel_od - 2*wall"),

    # ── section: outer liner ─────────────────────────────────────────────────
    ("outer_liner_od",
     fmt_mm(res['outer_liner_od_cold_mm']), "mm",
     "Outer liner OD — COLD BUILD dimension (room-temperature machining target). "
     "Hot (operating) = " + fmt_mm(res['outer_liner_od_hot_mm'])),

    ("outer_liner_id",
     fmt_mm(res['outer_liner_id_mm']), "mm",
     "Outer liner ID = outer_liner_od - 2*wall_thickness"),

    ("outer_annulus_gap",
     fmt_mm(res['outer_annulus_gap_mm']), "mm",
     "Radial width of outer feed annulus (casing ID to outer liner OD)/2"),

    # ── section: inner liner ─────────────────────────────────────────────────
    ("inner_liner_od",
     fmt_mm(res['inner_liner_od_cold_mm']), "mm",
     "Inner liner OD — COLD BUILD dimension. "
     "Hot = " + fmt_mm(res['inner_liner_od_hot_mm'])),

    ("inner_liner_id",
     fmt_mm(res['inner_liner_id_mm']), "mm",
     "Inner liner ID = inner_liner_od - 2*wall_thickness"),

    ("inner_liner_annulus_gap",
     fmt_mm(res['inner_annulus_gap_mm']), "mm",
     "Radial width of inner feed annulus (inner liner ID to shaft tunnel OD)/2"),

    # ── section: chamber length and zone lengths ──────────────────────────────
    ("chamber_length",
     fmt_mm(res['chamber_length_mm']), "mm",
     "Total axial combustion liner length"),

    ("prim_sec_length",
     fmt_mm(res['L_primary_mm']), "mm",
     "Primary zone axial length (30% of chamber_length). "
     "Secondary = " + fmt_mm(res['L_secondary_mm']) + " mm  "
     "— update this expression if your NX model uses prim_sec as a combined zone length."),

    ("dil_length",
     fmt_mm(res['L_dilution_mm']), "mm",
     "Dilution zone axial length (40% of chamber_length)"),

    # ── section: combustion annulus ───────────────────────────────────────────
    ("combustion_annulus_radial_height",
     fmt_mm(res['combustion_gap_mm']), "mm",
     "Radial height of the combustion annulus H = (outer_liner_id - inner_liner_od)/2"),

    # ── section: primary hole diameters ──────────────────────────────────────
    ("outer_primHole_dia",
     fmt_mm(res['pri_out_mm'], 2), "mm",
     "Outer liner primary zone hole diameter"),

    ("inner_primHole_dia",
     fmt_mm(res['pri_in_mm'], 2), "mm",
     "Inner liner primary zone hole diameter"),

    # ── section: secondary hole diameters ────────────────────────────────────
    ("outer_secHole_dia",
     fmt_mm(res['sec_out_mm'], 2), "mm",
     "Outer liner secondary zone hole diameter"),

    ("inner_secHole_dia",
     fmt_mm(res['sec_in_mm'], 2), "mm",
     "Inner liner secondary zone hole diameter"),

    # ── section: dilution hole diameters ─────────────────────────────────────
    ("outer_dilHole_dia",
     fmt_mm(res['dil_out_mm'], 2), "mm",
     "Outer liner dilution zone hole diameter"),

    ("inner_dilHole_dia",
     fmt_mm(res['dil_in_mm'], 2), "mm",
     "Inner liner dilution zone hole diameter"),

    # ── section: hole counts ─────────────────────────────────────────────────
    # V21 always produces equal counts on outer and inner liners, so one
    # expression drives both liners in the parametric model.
    # Inner counts: pri_in_qty={pri_in}, sec_in_qty={sec_in}, dil_in_qty={dil_in}
    # (identical to outer counts — included here for traceability).
    ("primHole_qty",
     fmt_int(res['pri_out_qty']), "int",
     f"Primary holes per liner (outer = inner = {res['pri_out_qty']})"),

    ("secHole_qty",
     fmt_int(res['sec_out_qty']), "int",
     f"Secondary holes per liner (outer = inner = {res['sec_out_qty']})"),

    ("dilHole_qty",
     fmt_int(res['dil_out_qty']), "int",
     f"Dilution holes per liner (outer = inner = {res['dil_out_qty']})"),

    # ── section: film cooling ─────────────────────────────────────────────────
    ("outer_film_hole_qty",
     fmt_int(res['film_holes_per_row_outer']), "int",
     f"Outer liner film cooling holes PER ROW ({res['film_n_rows']} rows × this = "
     f"{res['film_n_rows'] * res['film_holes_per_row_outer']} total outer film holes)"),

    ("inner_film_hole_qty",
     fmt_int(res['film_holes_per_row_inner']), "int",
     f"Inner liner film cooling holes PER ROW ({res['film_n_rows']} rows × this = "
     f"{res['film_n_rows'] * res['film_holes_per_row_inner']} total inner film holes)"),

    ("film_hole_area",
     fmt_mm2(res['film_total_area_mm2']), "mm2",
     f"Total film cooling hole area both liners combined "
     f"(hole dia {res['film_hole_dia_mm_actual']} mm, {res['film_total_holes']} holes total)"),
]

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Generate the NX journal
# ─────────────────────────────────────────────────────────────────────────────

timestamp  = datetime.now().strftime("%Y-%m-%d %H:%M")
out_name   = "Combustion_Journal_V2.py"
out_path   = os.path.join(SCRIPT_DIR, out_name)

# Build the expression block: FindObject → EditExpressionWithUnits for every entry
def build_expr_block(expressions):
    """
    Returns (variable_declarations, edit_calls, objects_list, n_exprs)
    so the journal can be written in the two-phase NX style:
      phase 1: all FindObject + EditExpression calls
      phase 2: MakeUpToDate(all_objects, markId) + DoUpdate
    """
    finds  = []   # expression_N = workPart.Expressions.FindObject("name")
    edits  = []   # workPart.Expressions.EditExpressionWithUnits(...)
    refs   = []   # expression_N  — for the objects array

    for i, (name, value, unit_tag, comment) in enumerate(expressions):
        var = f"expression{i+1}"
        refs.append(var)
        finds.append((var, name, comment))

        if unit_tag == 'mm':
            unit_ref = "unit_mm"
        elif unit_tag == 'mm2':
            unit_ref = "unit_mm2"
        else:
            unit_ref = "NXOpen.Unit.Null"

        edits.append((var, unit_ref, value, name))

    return finds, edits, refs

finds, edits, refs = build_expr_block(expressions)
n = len(expressions)

lines = []
add = lines.append

add('# NX Open Journal — Generated by GenerateNXJournal.py')
add(f'# Preset      : {preset_label}')
add(f'# Generated   : {timestamp}')
add(f'# Source tool : V21_CombustionChamberDesign.py')
add(f'# Key outputs : mdot_fuel={res["mdot_fuel"]*1000:.2f} g/s  '
    f'phi={res["overall_phi"]:.3f}  '
    f'gap={res["combustion_gap_mm"]:.1f}mm  '
    f'n_vap={res["vap_n"]}')
add('#')
add('# Inputs used:')
for k, v in chosen_inputs.items():
    add(f'#   {k} = {v}')
add('#')
add('# IMPORTANT — liner OD dimensions are COLD BUILD (room-temperature).')
add('# Parts will expand to their nominal hot dimensions at operating temperature.')
add('# Do NOT substitute the hot dimensions in NX — the thermal offset is already applied.')
add('')
add('import NXOpen')
add('')
add('')
add('def main():')
add('')
add('    theSession  = NXOpen.Session.GetSession()')
add('    workPart    = theSession.Parts.Work')
add('')
add('    # --- Unit references (fetched once, reused) ---------------------------')
add('    unit_mm  = workPart.UnitCollection.FindObject("MilliMeter")')
add('    unit_mm2 = workPart.UnitCollection.FindObject("SquareMilliMeter")')
add('')
add('    # --- Undo mark ---------------------------------------------------------')
add('    markId = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")')
add('    theSession.SetUndoMarkName(markId, "Combustor Expressions Update")')
add('')
add('    # --- Locate all expressions --------------------------------------------')

for var, name, comment in finds:
    add(f'    # {comment}')
    add(f'    {var} = workPart.Expressions.FindObject("{name}")')

add('')
add('    # --- Edit all expressions ----------------------------------------------')
add('    # Values sourced from V21_CombustionChamberDesign.py output.')

for var, unit_ref, value, name in edits:
    add(f'    workPart.Expressions.EditExpressionWithUnits({var}, {unit_ref}, "{value}")'
        f'  # {name}')

add('')
add('    # --- Commit: mark all for update, then regenerate the model -----------')
add(f'    objects = [NXOpen.NXObject.Null] * {n}')
for i, ref in enumerate(refs):
    add(f'    objects[{i}] = {ref}')
add('')
add('    markId_upd = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible,')
add('                                        "Make Up to Date")')
add('    theSession.UpdateManager.MakeUpToDate(objects, markId_upd)')
add('')
add('    markId_nx = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible,')
add('                                       "NX update")')
add('    nErrs = theSession.UpdateManager.DoUpdate(markId_nx)')
add('')
add('    theSession.DeleteUndoMark(markId_nx,  "NX update")')
add('    theSession.DeleteUndoMark(markId_upd, None)')
add('    theSession.SetUndoMarkName(markId, "Combustor Expressions Update")')
add('')
add('    if nErrs > 0:')
add('        print(f"  WARNING: NX reported {nErrs} update error(s). "')
add('              "Check model for failed features.")')
add('    else:')
add('        print("  Expressions updated and model regenerated successfully.")')
add('')
add('    # --- Save --------------------------------------------------------------')
add('    status = workPart.Save(NXOpen.BasePart.SaveComponents.TrueValue,')
add('                           NXOpen.BasePart.CloseAfterSave.FalseValue)')
add('    status.Dispose()')
add('    print("  Part saved.")')
add('')
add('')
add('if __name__ == "__main__":')
add('    main()')
add('')

journal_text = '\n'.join(lines)

with open(out_path, 'w', encoding='utf-8') as f:
    f.write(journal_text)

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Summary
# ─────────────────────────────────────────────────────────────────────────────

print()
print("=" * 65)
print(f"  Journal written to:")
print(f"    {out_path}")
print()
print(f"  {n} expressions will be updated:")
print()
col_w = 36
print(f"  {'NX expression':<{col_w}} {'New value':>12}  Unit")
print("  " + "-" * (col_w + 16))
for name, value, unit_tag, _ in expressions:
    unit_label = {"mm": "mm", "mm2": "mm²", "int": "integer"}[unit_tag]
    print(f"  {name:<{col_w}} {value:>12}  {unit_label}")

print()
print("  To use in NX:")
print("    Tools → Journals → Run…  → select Combustion_Journal_V2.py")
print("=" * 65)
print()
