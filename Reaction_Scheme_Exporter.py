# -*- coding: utf-8 -*-
"""
Reaction Scheme & Molecule Image Exporter
Composición del esquema con PIL a partir de cada molécula renderizada con RDKit.

Todas las magnitudes se definen en unidades físicas (mm, pt, cm, dpi), de modo que la imagen
exportada tiene un tamaño real conocido (p. ej. 8 cm = 1 columna, 16 cm = 2 columnas) y los
metadatos de resolución del PNG son correctos al insertarla en Word, LaTeX o Illustrator.
"""

import io
import re
import math
import hashlib
import streamlit as st

# --- Dependency Verification ---
try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor, rdChemReactions
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    from streamlit_ketcher import st_ketcher
    KETCHER_AVAILABLE = True
except ImportError:
    KETCHER_AVAILABLE = False

try:
    from PIL import Image, ImageDraw, ImageFont
    PIL_AVAILABLE = True
    RESAMPLE = getattr(Image, "Resampling", Image).LANCZOS
except ImportError:
    PIL_AVAILABLE = False


# ----------------------------------------------------------------------------
# Utilidades
# ----------------------------------------------------------------------------
def mm_to_px(mm, dpi):
    return mm * dpi / 25.4


def pt_to_px(pt, dpi):
    return pt * dpi / 72.0


def cm_to_px(cm, dpi):
    return int(round(cm / 2.54 * dpi))


def get_font(size, bold=False):
    """Devuelve una fuente TrueType disponible (Windows/Linux/macOS) o la de por defecto."""
    size = max(6, int(round(size)))
    candidates = (
        ["arialbd.ttf", "Arial Bold.ttf", "DejaVuSans-Bold.ttf", "Helvetica-Bold.ttf"]
        if bold
        else ["arial.ttf", "Arial.ttf", "DejaVuSans.ttf", "Helvetica.ttf"]
    )
    for name in candidates:
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            continue
    try:
        return ImageFont.load_default(size=size)  # Pillow >= 10.1
    except Exception:
        return ImageFont.load_default()


def text_size(draw, text, font):
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def draw_text_centered(draw, cx, top_y, text, font, fill=(0, 0, 0, 255)):
    """Dibuja texto centrado horizontalmente en cx, con su borde superior visible en top_y."""
    bbox = draw.textbbox((0, 0), text, font=font)
    w = bbox[2] - bbox[0]
    draw.text((cx - w / 2 - bbox[0], top_y - bbox[1]), text, fill=fill, font=font)


def wrap_text(draw, text, font, max_w):
    """
    Parte `text` en varias líneas (por espacios) para que no superen max_w px.
    Un espacio duro (\u00a0) impide el salto de línea en ese punto; se convierte en
    espacio normal en el resultado.
    """
    words = text.split(" ")
    words = [w for w in words if w]
    if not words:
        return []

    def width(t):
        return text_size(draw, t.replace("\u00a0", " "), font)[0]

    lines, cur = [], words[0]
    for w in words[1:]:
        trial = cur + " " + w
        if width(trial) <= max_w:
            cur = trial
        else:
            lines.append(cur)
            cur = w
    lines.append(cur)
    return [ln.replace("\u00a0", " ") for ln in lines]


def clean_condition(text):
    if text is None:
        return ""
    text = text.strip()
    return "" if text.upper() in ("", "N/A", "NA", "-") else text


def split_reaction(smiles_str):
    """
    Separa el SMILES de Ketcher en (hay_flecha, reactivos, agentes, productos).
    Cada elemento es una lista de SMILES individuales.
    """
    raw = (smiles_str or "").strip()
    smi = raw.split()[0] if raw else ""
    if not smi:
        return False, [], [], []

    if ">" in smi:
        parts = smi.split(">")
        reactants_part = parts[0]
        agents_part = parts[1] if len(parts) == 3 else ""
        products_part = parts[-1]
        has_arrow = True
    else:
        reactants_part, agents_part, products_part = smi, "", ""
        has_arrow = False

    def split_dots(x):
        return [t for t in x.split(".") if t.strip()]

    return has_arrow, split_dots(reactants_part), split_dots(agents_part), split_dots(products_part)



# ----------------------------------------------------------------------------
# Grupos genéricos (R1, R2…) y lectura de la salida de Ketcher
# ----------------------------------------------------------------------------
R_STYLES = {"R¹ (superíndice)": "sup", "R₁ (subíndice)": "sub", "R1 (texto plano)": "plain"}


def format_r_label(text, style="sub"):
    """'R1' -> 'R<sub>1</sub>' / 'R<sup>1</sup>' / 'R1'. Otros textos (Ar, X, Ph…) no se tocan."""
    text = (text or "").strip()
    m = re.fullmatch(r"([A-Za-z]+)(\d+)([a-z]?)", text)
    if not m or style == "plain":
        return text
    base, num, suf = m.groups()
    tag = "sup" if style == "sup" else "sub"
    return f"{base}<{tag}>{num}{suf}</{tag}>"


def count_generic_atoms(smiles_str):
    """(n_sin_etiqueta, n_con_etiqueta) de átomos genéricos (*) de una molécula en SMILES.
    Se consideran con etiqueta los que llevan número de mapeo, p. ej. [*:1] -> R1."""
    mol = Chem.MolFromSmiles((smiles_str or "").strip())
    if mol is None:
        mol = Chem.MolFromSmiles((smiles_str or "").strip(), sanitize=False)
    if mol is None:
        return 0, 0
    free = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and not a.GetAtomMapNum())
    lab = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetAtomMapNum())
    return free, lab


def apply_generic_labels(mol, smiles_str, rlabels=None, r_style="sub"):
    """
    Pone nombre a los átomos genéricos (*) para que RDKit no los dibuje como '*':
      - [*:n]  -> 'Rn'  (así se codifican los grupos R que vienen de un molfile)
      - '*' sin número -> texto que el usuario asignó en rlabels[smiles] (por orden de aparición)
    """
    overrides = (rlabels or {}).get((smiles_str or "").strip()) or []
    k = 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 0:
            continue
        text = None
        if a.GetAtomMapNum():
            text = f"R{a.GetAtomMapNum()}"
        else:
            if k < len(overrides) and str(overrides[k]).strip():
                text = str(overrides[k]).strip()
            k += 1
        if text:
            a.SetProp("atomLabel", format_r_label(text, r_style))
            a.SetAtomMapNum(0)
            a.SetIsotope(0)


def _mol_from_block(blk):
    mol = Chem.MolFromMolBlock(blk, sanitize=True, removeHs=True)
    if mol is None:
        mol = Chem.MolFromMolBlock(blk, sanitize=False, removeHs=False)
        if mol is not None:
            try:
                mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass
    return mol


def _frag_to_smiles(frag):
    """SMILES de un fragmento donde cada grupo R# del molfile queda codificado como [*:n]."""
    frag = Chem.Mol(frag)
    for a in frag.GetAtoms():
        if a.GetAtomicNum() == 0:
            n = a.GetIntProp("_MolFileRLabel") if a.HasProp("_MolFileRLabel") else 0
            a.SetIsotope(0)
            a.SetAtomMapNum(int(n) if n else 0)
    try:
        return Chem.MolToSmiles(frag)
    except Exception:
        return Chem.MolToSmiles(frag, canonical=False)


def _mol_to_smiles_list(mol):
    if mol is None or mol.GetNumAtoms() == 0:
        return []
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    return [_frag_to_smiles(f) for f in frags]


def parse_ketcher_output(text):
    """
    Convierte lo que devuelve Ketcher (SMILES, Molfile o RXN) en una cadena SMILES de reacción/
    moléculas que usa el resto de la aplicación. Los grupos R de un molfile/RXN se codifican como
    [*:n] y luego se dibujan como 'Rn'.
    Devuelve (smiles, info) con info = {"format": "smiles"|"molfile"|"rxn"|"error", "fragments": n}.
    """
    t = (text or "").strip()
    if not t:
        return "", {"format": "smiles", "fragments": 0}

    first = t.split("\n", 1)[0].strip()

    # ---- RXN ----
    if first.startswith("$RXN"):
        try:
            norm = t.replace("\r\n", "\n")
            lines = norm.split("\n")
            blocks = re.split(r"\$MOL\n", norm)[1:]
            reactants, products, agents = [], [], []
            ok = "V3000" not in first
            if ok:
                counts = lines[4].split()
                nr, npd = int(counts[0]), int(counts[1])
                na = int(counts[2]) if len(counts) > 2 else 0
                ok = len(blocks) == nr + npd + na
            if ok:
                mols = [_mol_from_block(b) for b in blocks]
                for i, m in enumerate(mols):
                    smi = _mol_to_smiles_list(m)
                    (reactants if i < nr else products if i < nr + npd else agents).extend(smi)
            else:  # RXN V3000 u otras variantes: lo lee RDKit
                rxn = rdChemReactions.ReactionFromRxnBlock(norm, sanitize=False)
                for m in rxn.GetReactants():
                    reactants.extend(_mol_to_smiles_list(m))
                for m in rxn.GetProducts():
                    products.extend(_mol_to_smiles_list(m))
                for m in rxn.GetAgents():
                    agents.extend(_mol_to_smiles_list(m))
            smi = f"{'.'.join(reactants)}>{'.'.join(agents)}>{'.'.join(products)}"
            return smi, {"format": "rxn", "fragments": len(reactants) + len(products) + len(agents)}
        except Exception:
            return "", {"format": "error", "fragments": 0}


    # ---- SMILES ----
    return t, {"format": "smiles", "fragments": 0}


def render_mol(smiles_str, bond_px=40.0, line_px=2, font_px=None, rlabels=None, r_style="sub"):
    """
    Renderiza una molécula con fondo transparente, con la longitud de enlace EXACTA en píxeles
    (bond_px), el grosor de línea line_px y el tamaño de letra de átomos font_px (ambos en px).
    Se recorta al bounding box real. Devuelve una imagen PIL RGBA o None si falla.
    """
    smiles_str = (smiles_str or "").strip()
    if not smiles_str:
        return None

    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles_str, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass

    try:
        apply_generic_labels(mol, smiles_str, rlabels, r_style)
    except Exception:
        pass

    try:
        rdDepictor.Compute2DCoords(mol)
        conf = mol.GetConformer()
        n = mol.GetNumAtoms()
        pts = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(n)]

        lens = [
            math.dist(pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()]) for b in mol.GetBonds()
        ]
        avg_len = (sum(lens) / len(lens)) if lens else 1.5
        k = bond_px / avg_len  # píxeles por unidad de coordenada

        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        margin_units = 2.5
        W = int((max(xs) - min(xs) + 2 * margin_units) * k) + 20
        H = int((max(ys) - min(ys) + 2 * margin_units) * k) + 20

        drawer = rdMolDraw2D.MolDraw2DCairo(W, H)
        opts = drawer.drawOptions()
        opts.fixedBondLength = float(k)
        opts.clearBackground = False
        opts.addAtomIndices = False
        opts.padding = 0.02
        opts.bondLineWidth = max(1, int(round(line_px)))
        if font_px:
            try:
                opts.fixedFontSize = int(round(font_px))
            except Exception:
                pass
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        img = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGBA")
        bbox = img.split()[3].getbbox()
        if not bbox:
            return None
        img = img.crop(bbox)
        if img.width < 3 or img.height < 3:
            return None
        return img
    except Exception:
        return None


# ----------------------------------------------------------------------------
# Composición del esquema
# ----------------------------------------------------------------------------
def build_scheme_image(
    smiles_str,
    solvent="",
    catalyst="",
    temperature="",
    labels=None,
    add_labels=True,
    bond_mm=5.0,
    dpi=600,
    line_pt=1.0,
    font_pt=8.0,
    number_style="(1)",
    start_number=1,
    rlabels=None,
    r_style="sub",
):
    has_arrow, reactants_smi, agents_smi, products_smi = split_reaction(smiles_str)
    if not (reactants_smi or products_smi):
        return None, [], []

    bl = mm_to_px(bond_mm, dpi)
    lw = max(1, int(round(pt_to_px(line_pt, dpi))))
    font_px = pt_to_px(font_pt, dpi)

    skipped = []

    reactants, products = [], []
    gidx = 0
    for sm in reactants_smi:
        img = render_mol(sm, bl, lw, font_px, rlabels, r_style)
        if img is None:
            skipped.append(sm)
        else:
            reactants.append((gidx, sm, img))
        gidx += 1
    for sm in products_smi:
        img = render_mol(sm, bl, lw, font_px, rlabels, r_style)
        if img is None:
            skipped.append(sm)
        else:
            products.append((gidx, sm, img))
        gidx += 1

    agents = []
    for sm in agents_smi:
        img = render_mol(sm, bl * 0.6, lw, font_px * 0.85, rlabels, r_style)
        if img is None:
            skipped.append(sm)
        else:
            agents.append((sm, img))

    if not reactants and not products:
        return None, [], skipped

    label_for = {}
    auto_n = int(start_number)
    for g, sm, img in reactants + products:
        if not add_labels:
            lbl = ""
        elif labels is None:
            lbl = f"({auto_n})" if number_style == "(1)" else f"{auto_n}"
            auto_n += 1
        else:
            lbl = (labels[g] if g < len(labels) else "").strip()
        label_for[g] = lbl

    PAD = int(round(0.6 * bl))
    ARROW_MARGIN = int(round(0.45 * bl))
    GAP_NUM = int(round(0.45 * bl))
    TXT_GAP = int(round(0.22 * bl))
    AGENT_GAP = int(round(0.3 * bl))
    LABEL_PAD = int(round(0.5 * bl))
    head = int(round(0.36 * bl))

    font_num = get_font(font_px, bold=True)
    font_plus = get_font(font_px * 1.35)
    font_cond = get_font(font_px * 0.95)

    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    plus_glyph_w = text_size(measure, "+", font_plus)[0]
    PLUS_W = plus_glyph_w + int(round(0.9 * bl))

    wrap_w = int(round(5.0 * bl))

    def cond_lines(texts):
        out = []
        for t in texts:
            t = clean_condition(t)
            if t:
                t = re.sub(r"(\d) (?=[A-Za-z°%µμ])", "\\1\u00a0", t)
                out.extend(wrap_text(measure, t, font_cond, wrap_w))
        return out

    above_lines = cond_lines([catalyst])
    below_lines = cond_lines([solvent, temperature])

    lh = text_size(measure, "Ag", font_cond)[1] + int(round(0.3 * font_px))

    agents_row_w = (
        sum(img.width for _, img in agents) + AGENT_GAP * (len(agents) - 1) if agents else 0
    )
    agents_row_h = max((img.height for _, img in agents), default=0)

    cond_widths = [text_size(measure, t, font_cond)[0] for t in above_lines + below_lines]
    arrow_w = max(
        int(round(2.8 * bl)),
        max(cond_widths, default=0) + int(round(0.6 * bl)),
        agents_row_w + int(round(0.6 * bl)),
    )

    placements = []
    label_items = []
    plus_centers = []

    def place_group(group, x):
        for i, (g, sm, img) in enumerate(group):
            lbl = label_for[g]
            lbl_w = text_size(measure, lbl, font_num)[0] + LABEL_PAD if lbl else 0
            slot = max(img.width, lbl_w)
            placements.append((img, x + (slot - img.width) / 2))
            label_items.append((x + slot / 2, lbl, sm))
            x += slot
            if i < len(group) - 1:
                plus_centers.append(x + PLUS_W / 2)
                x += PLUS_W
        return x

    x = place_group(reactants, PAD)

    arrow_x0 = arrow_x1 = None
    if has_arrow:
        x += ARROW_MARGIN
        arrow_x0 = x
        x += arrow_w
        arrow_x1 = x
        x += ARROW_MARGIN

    x = place_group(products, x)
    total_w = int(round(x + PAD))

    mol_h_max = max(img.height for _, _, img in reactants + products)
    half = mol_h_max / 2

    above_h = 0
    if has_arrow:
        above_h = lw + TXT_GAP + len(above_lines) * lh
        if agents:
            above_h += agents_row_h + AGENT_GAP
    below_h = (lw + TXT_GAP + len(below_lines) * lh) if has_arrow and below_lines else 0

    has_labels = any(label_for.values())
    num_h = text_size(measure, "Ay(1)|", font_num)[1] if has_labels else 0
    num_block = (half + GAP_NUM + num_h) if has_labels else half

    up = max(half, above_h)
    down = max(num_block, below_h, half)

    axis_y = PAD + up
    total_h = int(round(axis_y + down + PAD))

    canvas = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    for img, px in placements:
        py = int(round(axis_y - img.height / 2))
        canvas.paste(img, (int(round(px)), py), img)

    plus_h = text_size(draw, "+", font_plus)[1]
    for cx in plus_centers:
        draw_text_centered(draw, cx, axis_y - plus_h / 2, "+", font_plus)

    if has_arrow:
        draw.line([(arrow_x0, axis_y), (arrow_x1 - head, axis_y)], fill=(0, 0, 0, 255), width=lw)
        draw.polygon(
            [
                (arrow_x1, axis_y),
                (arrow_x1 - head, axis_y - head / 2.2),
                (arrow_x1 - head, axis_y + head / 2.2),
            ],
            fill=(0, 0, 0, 255),
        )
        arrow_cx = (arrow_x0 + arrow_x1) / 2

        cur = axis_y - lw - TXT_GAP
        for line in reversed(above_lines):
            cur -= lh
            draw_text_centered(draw, arrow_cx, cur, line, font_cond)
        if agents:
            cur -= AGENT_GAP
            ax = arrow_cx - agents_row_w / 2
            for _, img in agents:
                canvas.paste(img, (int(round(ax)), int(round(cur - img.height))), img)
                ax += img.width + AGENT_GAP

        cur = axis_y + lw + TXT_GAP
        for line in below_lines:
            draw_text_centered(draw, arrow_cx, cur, line, font_cond)
            cur += lh

    labeled = []
    if has_labels:
        num_y = axis_y + half + GAP_NUM
        for cx, lbl, sm in label_items:
            if lbl:
                draw_text_centered(draw, cx, num_y, lbl, font_num)
                labeled.append((lbl, sm))

    return canvas.convert("RGB"), labeled, skipped


def build_figure_image(
    smiles_str,
    labels=None,
    add_labels=True,
    bond_mm=5.0,
    dpi=600,
    line_pt=0.8,
    font_pt=8.0,
    number_style="(1)",
    start_number=1,
    cols=None,
    target_px=None,
    rlabels=None,
    r_style="sub",
):
    _, r_smi, a_smi, p_smi = split_reaction(smiles_str)
    mols_smi = r_smi + a_smi + p_smi
    if not mols_smi:
        return None, [], [], {}

    bl = mm_to_px(bond_mm, dpi)
    lw = max(1, int(round(pt_to_px(line_pt, dpi))))
    font_px = pt_to_px(font_pt, dpi)

    skipped = []
    items = []
    for g, sm in enumerate(mols_smi):
        img = render_mol(sm, bl, lw, font_px, rlabels, r_style)
        if img is None:
            skipped.append(sm)
        else:
            items.append((g, sm, img))
    if not items:
        return None, [], skipped, {}

    label_for = {}
    auto_n = int(start_number)
    for g, sm, img in items:
        if not add_labels:
            lbl = ""
        elif labels is None:
            lbl = f"({auto_n})" if number_style == "(1)" else f"{auto_n}"
            auto_n += 1
        else:
            lbl = (labels[g] if g < len(labels) else "").strip()
        label_for[g] = lbl

    PAD = int(round(0.6 * bl))
    GAP_X = int(round(0.9 * bl))
    GAP_Y = int(round(0.9 * bl))
    GAP_NUM = int(round(0.45 * bl))
    LABEL_PAD = int(round(0.5 * bl))

    font_num = get_font(font_px, bold=True)
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))

    slots = []
    for g, sm, img in items:
        lbl = label_for[g]
        lbl_w = text_size(measure, lbl, font_num)[0] + LABEL_PAD if lbl else 0
        slots.append(max(img.width, lbl_w))

    n = len(items)

    def col_widths(c):
        return [max(slots[j::c]) for j in range(c)]

    def grid_width(c):
        return 2 * PAD + sum(col_widths(c)) + (c - 1) * GAP_X

    if cols is None:
        cols = n
        if target_px:
            cols = 1
            for c in range(n, 0, -1):
                if grid_width(c) <= target_px:
                    cols = c
                    break
    cols = max(1, min(int(cols), n))
    rows = int(math.ceil(n / cols))
    cws = col_widths(cols)
    col_x = []
    xx = PAD
    for w_ in cws:
        col_x.append(xx)
        xx += w_ + GAP_X

    num_h_full = text_size(measure, "Ay(1)|", font_num)[1]

    row_groups = [items[i * cols:(i + 1) * cols] for i in range(rows)]
    row_mol_h, row_has_lbl, row_h = [], [], []
    for grp in row_groups:
        mh = max(img.height for _, _, img in grp)
        hl = any(label_for[g] for g, _, _ in grp)
        row_mol_h.append(mh)
        row_has_lbl.append(hl)
        row_h.append(mh + (GAP_NUM + num_h_full if hl else 0))

    total_w = int(grid_width(cols))
    total_h = int(2 * PAD + sum(row_h) + (rows - 1) * GAP_Y)

    canvas = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    labeled = []
    y = PAD
    for ri, grp in enumerate(row_groups):
        axis_y = y + row_mol_h[ri] / 2
        for j, (g, sm, img) in enumerate(grp):
            cx = col_x[j] + cws[j] / 2
            canvas.paste(img, (int(round(cx - img.width / 2)), int(round(axis_y - img.height / 2))), img)
            lbl = label_for[g]
            if lbl:
                draw_text_centered(draw, cx, y + row_mol_h[ri] + GAP_NUM, lbl, font_num)
                labeled.append((lbl, sm))
        y += row_h[ri] + GAP_Y

    return canvas.convert("RGB"), labeled, skipped, {"cols": cols, "rows": rows, "n": n}


def make_scheme(
    smiles_str,
    width_cm=8.0,
    dpi=600,
    fit_width=True,
    bond_mm=5.0,
    mode="reaction",
    cols=None,
    **kwargs,
):
    target_px = cm_to_px(width_cm, dpi)

    if mode == "figure":
        for k in ("solvent", "catalyst", "temperature"):
            kwargs.pop(k, None)
        fig_cols = cols
        if fit_width and cols is None:
            ref = build_figure_image(
                smiles_str, bond_mm=5.0, dpi=dpi, cols=None, target_px=target_px, **kwargs
            )
            if ref[0] is None:
                return None
            fig_cols = ref[3]["cols"]

        def build(b):
            return build_figure_image(
                smiles_str, bond_mm=b, dpi=dpi, cols=fig_cols, target_px=target_px, **kwargs
            )
    else:
        def build(b):
            img_, lab_, sk_ = build_scheme_image(smiles_str, bond_mm=b, dpi=dpi, **kwargs)
            return img_, lab_, sk_, {}

    if fit_width:
        b0 = 5.0
        img0, labeled, skipped, extra0 = build(b0)
        if img0 is None:
            return None
        tried = [(img0, labeled, skipped, extra0, b0)]
        w0 = img0.width
        b1 = min(12.0, max(1.0, b0 * target_px / w0))
        tol = max(2, int(0.003 * target_px))
        for _ in range(7):
            img1, labeled1, skipped1, extra1 = build(b1)
            tried.append((img1, labeled1, skipped1, extra1, b1))
            w1 = img1.width
            if abs(w1 - target_px) <= tol or b1 in (1.0, 12.0):
                break
            slope = (w1 - w0) / (b1 - b0) if abs(b1 - b0) > 1e-9 else 0
            b2 = b1 + (target_px - w1) / slope if slope > 0 else b1 * target_px / w1
            b0, w0 = b1, w1
            b1 = min(12.0, max(1.0, b2))
        close = [t for t in tried if abs(t[0].width - target_px) <= tol]
        fitting = [t for t in tried if t[0].width <= target_px]
        if close:
            best = min(close, key=lambda t: abs(t[0].width - target_px))
        elif fitting:
            best = max(fitting, key=lambda t: t[0].width)
        else:
            best = min(tried, key=lambda t: t[0].width)
        img, labeled, skipped, extra, b_eff = best
    else:
        img, labeled, skipped, extra = build(bond_mm)
        if img is None:
            return None
        b_eff = bond_mm

    exceeded = False
    if img.width > target_px:
        if fit_width:
            ratio = target_px / img.width
            img = img.resize((target_px, max(1, int(round(img.height * ratio)))), RESAMPLE)
            b_eff *= ratio
        else:
            exceeded = True
    elif img.width < target_px:
        padded = Image.new("RGB", (target_px, img.height), (255, 255, 255))
        padded.paste(img, ((target_px - img.width) // 2, 0))
        img = padded

    return {
        "image": img,
        "labeled": labeled,
        "skipped": skipped,
        "bond_mm": b_eff,
        "dpi": dpi,
        "width_cm": img.width / dpi * 2.54,
        "height_cm": img.height / dpi * 2.54,
        "exceeded": exceeded,
        "layout": extra,
    }


def to_png_bytes(img, dpi):
    buf = io.BytesIO()
    img.save(buf, format="PNG", dpi=(dpi, dpi))
    return buf.getvalue()


# ----------------------------------------------------------------------------
# Interfaz Streamlit
# ----------------------------------------------------------------------------
def _thumbnail(smiles_str, rlabels=None, r_style="sub"):
    img = render_mol(smiles_str, bond_px=22, line_px=1, font_px=11, rlabels=rlabels, r_style=r_style)
    if img is None:
        return None
    bg = Image.new("RGBA", (img.width + 16, img.height + 16), (255, 255, 255, 255))
    bg.paste(img, (8, 8), img)
    return bg.convert("RGB")


def render_reaction_scheme_exporter():
    st.title("🎨 Generador y Exportador de Esquemas Químicos")
    st.markdown(
        "Dibuja una **reacción química** o una **molécula**, introduce las **condiciones de reacción** y descárgala "
        "en formato PNG con alta resolución para tus informes, publicaciones o TFGs."
    )

    if not RDKIT_AVAILABLE or not KETCHER_AVAILABLE or not PIL_AVAILABLE:
        st.warning("Se requieren `RDKit`, `streamlit-ketcher` y `Pillow` para ejecutar este módulo.")
        return

    st.markdown("### 🧭 Tipo de imagen")
    image_type = st.radio(
        "¿Qué quieres generar?",
        ["Esquema de reacción", "Figura de moléculas"],
        horizontal=True,
        key="exp_image_type",
        help="Esquema: reactivos + flecha con condiciones + productos. "
             "Figura: moléculas alineadas en filas, sin signos + ni flecha.",
    )
    is_figure = image_type.startswith("Figura")

    st.markdown("### ✍️ Editor Molecular / Reacciones")
    if is_figure:
        st.caption(
            "Dibuja una o varias moléculas. Se colocarán alineadas en filas ajustadas al ancho elegido. "
            "Si dibujas una reacción, todos sus componentes se tratarán como moléculas independientes."
        )

    raw_output = ""
    try:
        raw_output = st_ketcher("", height=500)
    except Exception:
        raw_output = ""

    reaction_smiles, parse_info = parse_ketcher_output(raw_output)

    if parse_info.get("format") == "error":
        st.error("No se pudo interpretar lo que ha devuelto Ketcher.")
        return

    if not reaction_smiles or not reaction_smiles.strip():
        st.info("Dibuja la estructura o reacción en el lienzo y pulsa **Apply** para previsualizar.")
        return

    st.markdown("---")
    catalyst, solvent, temperature = "", "", ""
    if not is_figure:
        st.markdown("### ⚙️ Condiciones de Reacción")

        with st.container(border=True):
            col_c1, col_c2, col_c3 = st.columns(3)
            with col_c1:
                catalyst = st.text_input(
                    "Reactivos / catalizador (sobre la flecha):", value="N/A", key="exp_catalyst"
                )
            with col_c2:
                solvent = st.text_input("Disolvente (bajo la flecha):", value="Etanol", key="exp_solvent")
            with col_c3:
                temperature = st.text_input(
                    "Temperatura / Tiempo (bajo la flecha):", value="Reflujo, 2 h", key="exp_temp"
                )

    # ---------------- Grupos genéricos (R1, R2…) ----------------
    _, r_all, a_all, p_all = split_reaction(reaction_smiles)
    if is_figure:
        components = [("Molécula", i + 1, sm) for i, sm in enumerate(r_all + a_all + p_all)]
    else:
        components = (
            [("Reactivo", i + 1, sm) for i, sm in enumerate(r_all)]
            + [("Agente", i + 1, sm) for i, sm in enumerate(a_all)]
            + [("Producto", i + 1, sm) for i, sm in enumerate(p_all)]
        )

    rxn_hash = hashlib.md5(reaction_smiles.strip().encode("utf-8")).hexdigest()[:8]
    rlabels = {}
    r_style = "sub"
    free_groups = []
    any_generic = False
    seen_smiles = set()
    for kind, k, sm in components:
        if sm in seen_smiles:
            continue
        seen_smiles.add(sm)
        n_free, n_lab = count_generic_atoms(sm)
        if n_free or n_lab:
            any_generic = True
        if n_free:
            free_groups.append((kind, k, sm, n_free))

    if any_generic:
        st.markdown("### 🔣 Grupos genéricos (R₁, R₂…)")
        with st.container(border=True):
            style_choice = st.radio(
                "Formato de los grupos R:", list(R_STYLES.keys()), horizontal=True, key="exp_r_style"
            )
            r_style = R_STYLES[style_choice]

            if free_groups:
                st.caption(
                    "Hay grupos genéricos que Ketcher ha exportado como «*» (el formato SMILES no puede guardar "
                    "su nombre). Escribe el nombre de cada uno: se numeran por orden de aparición en el SMILES y la "
                    "miniatura se actualiza al momento, así puedes cambiarlos hasta que coincidan con tu dibujo. "
                    "Pueden ser R1, R2… o cualquier texto (Ar, X, Ph…). Si dibujaste R1, R2… en Ketcher, "
                    "usa el formato Molfile en «Opciones avanzadas del editor» para que se conserven solos."
                )
                for kind, k, sm, n_free in free_groups:
                    col_img, col_txt = st.columns([1, 3])
                    vals = []
                    with col_txt:
                        st.markdown(f"**{kind} {k}** · {n_free} grupo(s) genérico(s)")
                        sm_hash = hashlib.md5(sm.encode("utf-8")).hexdigest()[:6]
                        in_cols = st.columns(min(n_free, 4))
                        for i in range(n_free):
                            with in_cols[i % 4]:
                                vals.append(
                                    st.text_input(
                                        f"Grupo {i + 1}",
                                        value=f"R{i + 1}",
                                        key=f"exp_rg_{rxn_hash}_{sm_hash}_{i}",
                                    )
                                )
                    rlabels[sm] = vals
                    with col_img:
                        thumb = _thumbnail(sm, rlabels, r_style)
                        if thumb is not None:
                            st.image(thumb, width=min(thumb.width, 150))
            else:
                st.caption("Grupos R detectados y etiquetados automáticamente (formato Molfile).")

    # ---------------- Etiquetas de los compuestos ----------------
    st.markdown("### 🔢 Etiquetas de las moléculas" if is_figure else "### 🔢 Etiquetas de los compuestos")

    _, reactants_smi, agents_smi, products_smi = split_reaction(reaction_smiles)
    if is_figure:
        all_smi = reactants_smi + agents_smi + products_smi
        compounds = [("Molécula", i + 1, sm) for i, sm in enumerate(all_smi)]
    else:
        compounds = [("Reactivo", i + 1, sm) for i, sm in enumerate(reactants_smi)] + [
            ("Producto", i + 1, sm) for i, sm in enumerate(products_smi)
        ]
    type_tag = "fig" if is_figure else "rxn"

    labels = None
    number_style, start_number = "(1)", 1

    with st.container(border=True):
        add_labels = st.checkbox(
            "Añadir etiqueta bajo cada molécula" if is_figure else "Añadir etiqueta bajo cada compuesto", value=True, key="exp_add_numbers"
        )

        if add_labels:
            mode = st.radio(
                "Modo de etiquetado:",
                ["Automático (1, 2, 3…)", "Personalizado"],
                horizontal=True,
                key="exp_label_mode",
            )

            if mode.startswith("Automático"):
                col_n1, col_n2 = st.columns(2)
                with col_n1:
                    number_style = st.selectbox("Formato del número:", ["(1)", "1"], key="exp_num_style")
                with col_n2:
                    start_number = int(
                        st.number_input("Empezar en:", min_value=1, max_value=999, value=1, step=1, key="exp_num_start")
                    )
            else:
                st.caption(
                    "Escribe el texto que quieras bajo cada compuesto (p. ej. `1a`, `(2)`, `Anilina`). "
                    "Déjalo vacío para no poner nada en ese compuesto."
                )
                labels = []
                for g, (kind, k, sm) in enumerate(compounds):
                    col_img, col_txt = st.columns([1, 3])
                    with col_img:
                        thumb = _thumbnail(sm, rlabels, r_style)
                        if thumb is not None:
                            st.image(thumb, width=min(thumb.width, 150))
                        else:
                            st.caption(f"`{sm}`")
                    with col_txt:
                        val = st.text_input(
                            f"{kind} {k}",
                            value=f"({g + 1})",
                            key=f"exp_lbl_{type_tag}_{rxn_hash}_{g}",
                        )
                    labels.append(val)

    # ---------------- Tamaño y exportación ----------------
    st.markdown("---")
    st.markdown("### 📐 Tamaño de la imagen")

    with st.container(border=True):
        col_w, col_cm, col_dpi = st.columns([2, 1, 1])
        with col_w:
            width_choice = st.selectbox(
                "Ancho de la imagen:",
                ["8 cm (1 columna)", "16 cm (2 columnas)", "Personalizado"],
                key="exp_width_choice",
            )
        with col_cm:
            if width_choice.startswith("8"):
                width_cm = 8.0
            elif width_choice.startswith("16"):
                width_cm = 16.0
            else:
                width_cm = float(
                    st.number_input(
                        "Ancho (cm):", min_value=3.0, max_value=40.0, value=12.0, step=0.5, key="exp_width_cm"
                    )
                )
            if width_choice != "Personalizado":
                st.markdown(f"**{width_cm:g} cm**")
        with col_dpi:
            dpi = int(st.selectbox("Resolución (dpi):", [300, 600, 1200], index=1, key="exp_dpi"))

        size_mode = st.radio(
            "Tamaño de las moléculas dentro de la imagen:" if is_figure else "Tamaño del esquema dentro de la imagen:",
            ["Ajustar al ancho (recomendado)", "Longitud de enlace fija"],
            horizontal=True,
            key="exp_size_mode",
        )
        fit_width = size_mode.startswith("Ajustar")

        cols = None
        if is_figure:
            cols_choice = st.selectbox(
                "Moléculas por fila:",
                ["Automático"] + list(range(1, 9)),
                key="exp_cols",
                help="Automático: tantas como quepan. Las moléculas se alinean en columnas "
                     "y pasan a la fila siguiente cuando no caben.",
            )
            cols = None if cols_choice == "Automático" else int(cols_choice)

        col_b, col_l, col_f = st.columns(3)
        with col_b:
            bond_mm = st.slider(
                "Longitud de enlace (mm):",
                min_value=2.0, max_value=15.0, value=5.0, step=0.5,
                key="exp_bond_mm", disabled=fit_width,
                help="Solo en modo «Longitud de enlace fija». ACS usa ≈ 5,1 mm (14,4 pt).",
            )
        with col_l:
            line_pt = st.slider(
                "Grosor de línea (pt):", min_value=0.5, max_value=3.0, value=0.8, step=0.1, key="exp_line_pt"
            )
        with col_f:
            font_pt = st.slider(
                "Tamaño de texto (pt):", min_value=6.0, max_value=24.0, value=8.0, step=0.5, key="exp_font_pt",
                help="Afecta a las letras de los átomos, las etiquetas y las condiciones de reacción."
            )

    st.markdown("### 🖼️ Previsualización y Descarga")

    default_name = "Figura_Moleculas" if is_figure else "Esquema_Reaccion"
    file_name_input = st.text_input(
        "Nombre del archivo (sin extensión):", value=default_name, key=f"exp_filename_{type_tag}"
    )

    save_editable_format = st.checkbox(
        "Mostrar botón para descargar la estructura química editable (formato Ketcher/SMILES)", 
        value=True, 
        key="exp_save_editable"
    )

    try:
        result = make_scheme(
            reaction_smiles,
            width_cm=width_cm,
            dpi=dpi,
            fit_width=fit_width,
            bond_mm=bond_mm,
            mode="figure" if is_figure else "reaction",
            cols=cols,
            solvent=solvent,
            catalyst=catalyst,
            temperature=temperature,
            labels=labels,
            add_labels=add_labels,
            line_pt=line_pt,
            font_pt=font_pt,
            number_style=number_style,
            start_number=start_number,
            rlabels=rlabels,
            r_style=r_style,
        )
    except Exception as e:
        result = None
        st.error(f"Error al generar la imagen: {e}")

    if result:
        if result["skipped"]:
            st.warning(
                "No se pudieron dibujar estas estructuras: "
                + ", ".join(f"`{x}`" for x in result["skipped"])
            )

        img_bytes = to_png_bytes(result["image"], result["dpi"])
        st.image(img_bytes, use_container_width=True)

        st.caption(
            f"Tamaño exportado: **{result['width_cm']:.1f} × {result['height_cm']:.1f} cm** "
            f"({result['image'].width} × {result['image'].height} px a {result['dpi']} dpi) · "
            f"longitud de enlace ≈ **{result['bond_mm']:.1f} mm**"
        )
        lay = result.get("layout") or {}
        if is_figure and lay:
            st.caption(
                f"{lay['n']} molécula(s) en {lay['rows']} fila(s) × {lay['cols']} columna(s)."
            )
        if result["exceeded"]:
            st.warning(
                f"Con esta longitud de enlace el contenido no cabe en {width_cm:g} cm; se exporta con "
                f"{result['width_cm']:.1f} cm de ancho. Reduce la longitud de enlace o usa «Ajustar al ancho»."
            )
        if fit_width and result["bond_mm"] < 3.0:
            st.warning(
                (
                    "Con este ancho y número de columnas los enlaces quedan muy cortos. "
                    "Prueba con menos moléculas por fila, un ancho mayor o etiquetas más cortas."
                    if is_figure
                    else "El esquema es muy largo para este ancho y los enlaces quedan muy cortos. "
                         "Considera 16 cm, menos compuestos o etiquetas más cortas."
                )
            )

        if result["labeled"]:
            with st.expander("📋 Correspondencia etiqueta ↔ compuesto (SMILES)"):
                for lbl, sm in result["labeled"]:
                    st.markdown(f"**{lbl}** — `{sm}`")

        clean_file_name = file_name_input.strip() or default_name
        full_png_name = (
            clean_file_name if clean_file_name.lower().endswith(".png") else f"{clean_file_name}.png"
        )
        full_smi_name = (
            clean_file_name if clean_file_name.lower().endswith(".smi") else f"{clean_file_name}.smi"
        )

        if save_editable_format:
            col_btn1, col_btn2 = st.columns(2)
            with col_btn1:
                st.download_button(
                    label="💾 Descargar Imagen PNG",
                    data=img_bytes,
                    file_name=full_png_name,
                    mime="image/png",
                    type="primary",
                    use_container_width=True,
                )
            with col_btn2:
                st.download_button(
                    label="🧪 Descargar Estructura Editable (.smi)",
                    data=reaction_smiles.encode("utf-8"),
                    file_name=full_smi_name,
                    mime="text/plain",
                    type="secondary",
                    use_container_width=True,
                    help="Guarda la cadena de la reacción/molécula. Puedes volver a importarla para editarla en el futuro."
                )
        else:
            st.download_button(
                label="💾 Descargar Imagen PNG",
                data=img_bytes,
                file_name=full_png_name,
                mime="image/png",
                type="primary",
                use_container_width=True,
            )
    elif result is None:
        st.error("No se pudo generar la imagen. Verifica la estructura dibujada.")

# Para uso standalone:
# if __name__ == "__main__":
    # render_reaction_scheme_exporter()