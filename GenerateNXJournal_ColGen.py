import os
import sys
import random
import importlib.util
from datetime import datetime
#
LINER_GROUP = [
    "EXTRUDE(31)",   # outer liner
    "EXTRUDE(6)",    # inner liner
    "EXTRUDE(4)",    # front plate
]

INDIVIDUAL_BODIES = {
    "casing":       "EXTRUDE(2)",    # outer casing
    "shaft_tunnel": "EXTRUDE(8)",    # shaft tunnel
}


NX_PALETTE = [
    # Reds / oranges
     11,   # vivid red
     16,   # brick red
     21,   # orange-red
     26,   # burnt orange
    # Yellows / golds
     31,   # golden yellow
     36,   # amber
     41,   # yellow-green
    # Greens
     46,   # olive
     51,   # mid green
     56,   # leaf green
     61,   # teal green
    # Cyans / aquas
     66,   # aqua
     71,   # cyan
     76,   # sky cyan
    # Blues
     81,   # sky blue
     86,   # cornflower blue
     91,   # medium blue
     96,   # royal blue
    101,   # navy blue
    # Purples / violets
    106,   # indigo
    111,   # purple
    116,   # violet
    # Pinks / magentas
    121,   # magenta
    126,   # hot pink
    131,   # rose pink
    # Muted / pastel (useful for backgrounds / less dominant parts)
    151,   # powder blue
    156,   # pale teal
    161,   # mint
    171,   # pale lime
    176,   # peach / salmon
    181,   # tan / warm beige
    # Grays (use sparingly — reserve for structural/background parts)
    191,   # medium gray
    196,   # silver-gray
    201,   # light silver
]

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
    sys.exit(1)

spec = importlib.util.spec_from_file_location("V21_combustor", v21_path)
v21  = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v21)

MicroJetCombustor = v21.MicroJetCombustor
user_inputs_6in   = v21.user_inputs_6in
kj66_inputs       = v21.kj66_inputs
KEYS              = v21.KEYS

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

try:
    model = MicroJetCombustor(chosen_inputs)
    res   = model.run()
except ValueError as e:
    print(f"\n  GEOMETRY ERROR: {e}")
    sys.exit(1)

print(f"\n  Model run complete.")
print(f"    mdot_fuel : {res['mdot_fuel']*1000:.2f} g/s")
print(f"    phi       : {res['overall_phi']:.3f}")
print(f"    gap       : {res['combustion_gap_mm']:.1f} mm")
print(f"    n_vap     : {res['vap_n']}")


available = NX_PALETTE.copy()
random.shuffle(available)

# Number of color slots needed: 1 for liner group + 1 per individual body
n_slots = 1 + len(INDIVIDUAL_BODIES)
if n_slots > len(available):
    # Wrap around if model has more parts than curated colors
    available = available * (n_slots // len(available) + 1)

liner_color = available.pop(0)

individual_colors = {}   # label → color int
for label in INDIVIDUAL_BODIES:
    individual_colors[label] = available.pop(0)

def fmt_mm(v, dp=3):  return f"{round(float(v), dp):.{dp}f}"
def fmt_int(v):       return str(int(round(float(v))))
def fmt_mm2(v, dp=1): return f"{round(float(v), dp):.{dp}f}"

expressions = [
    ("casing_od",                        fmt_mm(res['casing_od_mm']),              "mm",  "Outer casing OD"),
    ("casing_id",                        fmt_mm(res['casing_id_mm']),              "mm",  "Outer casing ID"),
    ("wall_thickness",                   fmt_mm(chosen_inputs['wall_thickness_mm']),"mm", "Nominal wall thickness"),
    ("shaft_tunnel_od",                  fmt_mm(res['shaft_tunnel_od_mm']),         "mm", "Shaft tunnel OD"),
    ("shaft_tunnel_id",                  fmt_mm(res['shaft_tunnel_id_mm']),         "mm", "Shaft tunnel ID"),
    ("outer_liner_od",                   fmt_mm(res['outer_liner_od_cold_mm']),     "mm",
     f"Outer liner OD — COLD BUILD. Hot = {fmt_mm(res['outer_liner_od_hot_mm'])}"),
    ("outer_liner_id",                   fmt_mm(res['outer_liner_id_mm']),          "mm", "Outer liner ID"),
    ("outer_annulus_gap",                fmt_mm(res['outer_annulus_gap_mm']),       "mm", "Outer feed annulus radial gap"),
    ("inner_liner_od",                   fmt_mm(res['inner_liner_od_cold_mm']),     "mm",
     f"Inner liner OD — COLD BUILD. Hot = {fmt_mm(res['inner_liner_od_hot_mm'])}"),
    ("inner_liner_id",                   fmt_mm(res['inner_liner_id_mm']),          "mm", "Inner liner ID"),
    ("inner_liner_annulus_gap",          fmt_mm(res['inner_annulus_gap_mm']),       "mm", "Inner feed annulus radial gap"),
    ("chamber_length",                   fmt_mm(res['chamber_length_mm']),          "mm", "Total liner axial length"),
    ("prim_sec_length",                  fmt_mm(res['L_primary_mm']),               "mm",
     f"Primary zone length (30%). Secondary = {fmt_mm(res['L_secondary_mm'])} mm"),
    ("dil_length",                       fmt_mm(res['L_dilution_mm']),              "mm", "Dilution zone length (40%)"),
    ("combustion_annulus_radial_height", fmt_mm(res['combustion_gap_mm']),          "mm", "Combustion annulus radial height H"),
    ("outer_primHole_dia",               fmt_mm(res['pri_out_mm'], 2),              "mm", "Outer liner primary hole dia"),
    ("inner_primHole_dia",               fmt_mm(res['pri_in_mm'],  2),              "mm", "Inner liner primary hole dia"),
    ("outer_secHole_dia",                fmt_mm(res['sec_out_mm'], 2),              "mm", "Outer liner secondary hole dia"),
    ("inner_secHole_dia",                fmt_mm(res['sec_in_mm'],  2),              "mm", "Inner liner secondary hole dia"),
    ("outer_dilHole_dia",                fmt_mm(res['dil_out_mm'], 2),              "mm", "Outer liner dilution hole dia"),
    ("inner_dilHole_dia",                fmt_mm(res['dil_in_mm'],  2),              "mm", "Inner liner dilution hole dia"),
    ("primHole_qty",    fmt_int(res['pri_out_qty']),  "int", f"Primary holes/liner ({res['pri_out_qty']} each)"),
    ("secHole_qty",     fmt_int(res['sec_out_qty']),  "int", f"Secondary holes/liner ({res['sec_out_qty']} each)"),
    ("dilHole_qty",     fmt_int(res['dil_out_qty']),  "int", f"Dilution holes/liner ({res['dil_out_qty']} each)"),
    ("outer_film_hole_qty", fmt_int(res['film_holes_per_row_outer']), "int",
     f"Outer film holes/row ({res['film_n_rows']} rows)"),
    ("inner_film_hole_qty", fmt_int(res['film_holes_per_row_inner']), "int",
     f"Inner film holes/row ({res['film_n_rows']} rows)"),
    ("film_hole_area",  fmt_mm2(res['film_total_area_mm2']), "mm2",
     f"Total film cooling area — {res['film_total_holes']} holes × {res['film_hole_dia_mm_actual']} mm dia"),
]

n_exprs = len(expressions)

timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
out_name  = "Combustion_Journal_V2.py"
out_path  = os.path.join(SCRIPT_DIR, out_name)

L = []
w = L.append

w('# NX Open Journal — Generated by GenerateNXJournalWithColors.py')
w(f'# Preset    : {preset_label}')
w(f'# Generated : {timestamp}')
w(f'# Sections  : [1] Expression updates  [2] Part color randomization')
w(f'#')
w(f'# Key outputs : mdot_fuel={res["mdot_fuel"]*1000:.2f} g/s  '
  f'phi={res["overall_phi"]:.3f}  '
  f'gap={res["combustion_gap_mm"]:.1f} mm  '
  f'n_vap={res["vap_n"]}')
w('#')
w('# Inputs used:')
for k, v in chosen_inputs.items():
    w(f'#   {k} = {v}')
w('#')
w('# Color assignments (rolled at generation time — re-run generator for new colors):')
w(f'#   Liner group  (outer_liner + inner_liner + front_liner_plates) → NX color {liner_color}')
for label, color in individual_colors.items():
    w(f'#   {label:<30} → NX color {color}')
w('#')
w('# LINER OD NOTE: outer_liner_od / inner_liner_od are COLD BUILD dimensions.')
w('# The model will expand to the nominal hot dimensions at operating temperature.')
w('')
w('import NXOpen')
w('')
w('')
w('def main():')
w('')
w('    theSession  = NXOpen.Session.GetSession()')
w('    workPart    = theSession.Parts.Work')
w('')

# ── Section 1: Expression updates ────────────────────────────────────────────
w('    # ══════════════════════════════════════════════════════════════════════')
w('    # SECTION 1 — EXPRESSION UPDATES')
w('    # ══════════════════════════════════════════════════════════════════════')
w('')
w('    unit_mm  = workPart.UnitCollection.FindObject("MilliMeter")')
w('    unit_mm2 = workPart.UnitCollection.FindObject("SquareMilliMeter")')
w('')
w('    markId_expr = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")')
w('    theSession.SetUndoMarkName(markId_expr, "Combustor Expressions Update")')
w('')
w('    # --- Locate all expressions ---')
for i, (name, value, unit_tag, comment) in enumerate(expressions):
    w(f'    # {comment}')
    w(f'    expression{i+1} = workPart.Expressions.FindObject("{name}")')
w('')
w('    # --- Edit all expressions ---')
for i, (name, value, unit_tag, comment) in enumerate(expressions):
    if unit_tag == 'mm':
        unit_ref = 'unit_mm'
    elif unit_tag == 'mm2':
        unit_ref = 'unit_mm2'
    else:
        unit_ref = 'NXOpen.Unit.Null'
    w(f'    workPart.Expressions.EditExpressionWithUnits(expression{i+1}, {unit_ref}, "{value}")  # {name}')
w('')
w(f'    objects_expr = [NXOpen.NXObject.Null] * {n_exprs}')
for i in range(n_exprs):
    w(f'    objects_expr[{i}] = expression{i+1}')
w('')
w('    markId_upd = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Make Up to Date")')
w('    theSession.UpdateManager.MakeUpToDate(objects_expr, markId_upd)')
w('    markId_nx  = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "NX update")')
w('    nErrs_expr = theSession.UpdateManager.DoUpdate(markId_nx)')
w('    theSession.DeleteUndoMark(markId_nx,  "NX update")')
w('    theSession.DeleteUndoMark(markId_upd, None)')
w('    theSession.SetUndoMarkName(markId_expr, "Combustor Expressions Update")')
w('')
w('    if nErrs_expr > 0:')
w('        print(f"  WARNING: {nErrs_expr} expression update error(s). Check model for failed features.")')
w('    else:')
w('        print("  [1/2] Expressions updated successfully.")')
w('')

# ── Section 2: Color assignments ─────────────────────────────────────────────
w('    # ══════════════════════════════════════════════════════════════════════')
w('    # SECTION 2 — COLOR ASSIGNMENTS')
w('    # ══════════════════════════════════════════════════════════════════════')
w('#')
w('#  If a FindObject call raises an exception, that body name does not exist')
w('#  in the current model.  Update LINER_GROUP / INDIVIDUAL_BODIES in the')
w('#  generator script (GenerateNXJournalWithColors.py) and re-generate.')
w('')

# Helper that writes one color-apply block
mark_counter = [100]   # mutable counter — avoids markId name collisions

def write_color_block(label, body_vars, color, comment=""):
    """
    Emits the FindObject + DisplayModification block for a list of body variable names.
    body_vars is a list of (var_name, nx_body_name) tuples.
    """
    mc = mark_counter[0]
    mark_counter[0] += 1
    mod_var = f'dispMod_{mc}'
    objects_var = f'objects_color_{mc}'
    mark_var = f'markId_color_{mc}'
    n_bodies = len(body_vars)

    # FindObject calls
    for var, nx_name in body_vars:
        w(f'    # {label} — body: {nx_name}')
        w(f'    {var} = workPart.Bodies.FindObject("{nx_name}")')

    w(f'    {mark_var} = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible,')
    w(f'                                        "Edit Object Display")')
    w(f'    {mod_var} = theSession.DisplayManager.NewDisplayModification()')
    w(f'    {mod_var}.ApplyToAllFaces    = True')
    w(f'    {mod_var}.ApplyToOwningParts = False')
    w(f'    {mod_var}.NewColor = {color}  # NX palette index')
    w(f'    {objects_var} = [NXOpen.DisplayableObject.Null] * {n_bodies}')
    for idx, (var, _) in enumerate(body_vars):
        w(f'    {objects_var}[{idx}] = {var}')
    w(f'    {mod_var}.Apply({objects_var})')
    w(f'    theSession.UpdateManager.DoUpdate({mark_var})')
    w(f'    {mod_var}.Dispose()')
    w('')

# Liner group — all bodies in one Apply() call (same color)
liner_body_vars = [(f'body_liner_{i}', nx_name)
                   for i, nx_name in enumerate(LINER_GROUP)]

w(f'    # ── Liner group: outer_liner + inner_liner + front_liner_plate(s)')
w(f'    # ── All receive the same color (NX index {liner_color}) ──────────────')
w('    try:')
# Indent the block
temp = []
saved_w = w.__self__ if hasattr(w, '__self__') else None

# We need indented code inside try/except — easier to build separately
liner_lines = []
def lw(s): liner_lines.append(s)

for var, nx_name in liner_body_vars:
    lw(f'        body_liner_{liner_body_vars.index((var, nx_name))} = workPart.Bodies.FindObject("{nx_name}")')

mc = mark_counter[0]; mark_counter[0] += 1
mod_var      = f'dispMod_{mc}'
objects_var  = f'objects_color_{mc}'
mark_var     = f'markId_color_{mc}'
n_liners     = len(LINER_GROUP)

lw(f'        {mark_var} = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Object Display")')
lw(f'        {mod_var} = theSession.DisplayManager.NewDisplayModification()')
lw(f'        {mod_var}.ApplyToAllFaces    = True')
lw(f'        {mod_var}.ApplyToOwningParts = False')
lw(f'        {mod_var}.NewColor = {liner_color}  # liner group color')
lw(f'        {objects_var} = [NXOpen.DisplayableObject.Null] * {n_liners}')
for idx, (var, _) in enumerate(liner_body_vars):
    lw(f'        {objects_var}[{idx}] = {var}')
lw(f'        {mod_var}.Apply({objects_var})')
lw(f'        theSession.UpdateManager.DoUpdate({mark_var})')
lw(f'        {mod_var}.Dispose()')
lw(f'        print("  Liner group color applied  (NX {liner_color})")')

for ll in liner_lines:
    w(ll)
w('    except NXOpen.NXException as e:')
w('        print(f"  WARNING — liner group: could not find one or more bodies: {e}")')
w('        print("  Check LINER_GROUP names in GenerateNXJournalWithColors.py")')
w('')

# Individual bodies — each in its own try/except
w('    # ── Individual bodies — each gets its own color ───────────────────────')
for label, nx_name in INDIVIDUAL_BODIES.items():
    color = individual_colors[label]
    mc = mark_counter[0]; mark_counter[0] += 1
    var         = f'body_{label}'
    mod_var     = f'dispMod_{mc}'
    objects_var = f'objects_color_{mc}'
    mark_var    = f'markId_color_{mc}'

    w(f'    # {label} → NX color {color}')
    w(f'    try:')
    w(f'        {var} = workPart.Bodies.FindObject("{nx_name}")')
    w(f'        {mark_var} = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Object Display")')
    w(f'        {mod_var} = theSession.DisplayManager.NewDisplayModification()')
    w(f'        {mod_var}.ApplyToAllFaces    = True')
    w(f'        {mod_var}.ApplyToOwningParts = False')
    w(f'        {mod_var}.NewColor = {color}')
    w(f'        objects_color_{mc} = [NXOpen.DisplayableObject.Null] * 1')
    w(f'        objects_color_{mc}[0] = {var}')
    w(f'        {mod_var}.Apply(objects_color_{mc})')
    w(f'        theSession.UpdateManager.DoUpdate({mark_var})')
    w(f'        {mod_var}.Dispose()')
    w(f'        print("  Color applied: {label}  (NX {color})")')
    w(f'    except NXOpen.NXException as e:')
    w(f'        print(f"  WARNING — {label}: body \\"{nx_name}\\" not found: {{e}}")')
    w('')

# ── Save ──────────────────────────────────────────────────────────────────────
w('    # ── Save ──────────────────────────────────────────────────────────────')
w('    status = workPart.Save(NXOpen.BasePart.SaveComponents.TrueValue,')
w('                           NXOpen.BasePart.CloseAfterSave.FalseValue)')
w('    status.Dispose()')
w('    print("  [2/2] Colors applied. Part saved.")')
w('')
w('')
w('if __name__ == "__main__":')
w('    main()')
w('')

journal_text = '\n'.join(L)

with open(out_path, 'w', encoding='utf-8') as f:
    f.write(journal_text)

print(f"  Journal written to:")
print(f"    {out_path}")
print()
print(f"  {n_exprs} expressions")

print(f"    {'Liner group (outer+inner+plates)':<40} --> NX {liner_color}")
for label, color in individual_colors.items():
    nx_name = INDIVIDUAL_BODIES[label]
    print(f"    {label:<40} → NX {color}")
print()
