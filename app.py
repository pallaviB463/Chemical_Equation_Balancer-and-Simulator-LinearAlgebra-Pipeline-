import streamlit as st
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sympy import Matrix, lcm
from math import gcd
from functools import reduce

st.set_page_config(page_title="Chemical Equation Balancer", layout="wide")
st.title("🧪 Advanced Chemical Equation Balancer & Simulator")

# ════════════════════════════════════════════════════════
# TIER 1 — FEATURE 1: RECURSIVE FORMULA PARSER
# Handles parentheses: Ca(OH)2, Al2(SO4)3, Fe2(SO4)3
# ════════════════════════════════════════════════════════

def parse_compound(formula):
    """
    Recursive descent parser for chemical formulas.
    Handles: simple (H2O), multi-letter elements (Fe, Ca),
             parentheses (Ca(OH)2), nested groups (Al2(SO4)3).
    Returns dict {element: count}.
    """
    def parse_tokens(s, i=0):
        counts = {}
        while i < len(s):
            if s[i] == '(':
                sub, i = parse_tokens(s, i + 1)
                num_start = i
                while i < len(s) and s[i].isdigit():
                    i += 1
                multiplier = int(s[num_start:i]) if i > num_start else 1
                for el, cnt in sub.items():
                    counts[el] = counts.get(el, 0) + cnt * multiplier
            elif s[i] == ')':
                return counts, i + 1
            elif s[i].isupper():
                j = i + 1
                while j < len(s) and s[j].islower():
                    j += 1
                el = s[i:j]
                i = j
                num_start = i
                while i < len(s) and s[i].isdigit():
                    i += 1
                num = int(s[num_start:i]) if i > num_start else 1
                counts[el] = counts.get(el, 0) + num
            else:
                i += 1
        return counts, i

    result, _ = parse_tokens(formula)
    return result


# ════════════════════════════════════════════════════════
# FEATURE 2: MINIMAL INTEGER COEFFICIENTS (GCD REDUCTION)
# e.g. [4,14,8,12] / GCD(4,14,8,12)=2 → [2,7,4,6]
# ════════════════════════════════════════════════════════

def solve(matrix):
    """
    1. Compute null space (SymPy exact arithmetic — no float errors)
    2. Clear fractions: multiply by LCM of all denominators
    3. Reduce to minimal: divide by GCD of all coefficients
    Returns (minimal_int_coeffs, raw_null_vector, full_nullspace)
    """
    ns = matrix.nullspace()
    if not ns:
        return None, None, ns
    sol = ns[0]
    lcm_val = lcm([term.q for term in sol])
    coeffs  = [int(term * lcm_val) for term in sol]
    g = reduce(gcd, [abs(c) for c in coeffs])
    coeffs = [c // g for c in coeffs]
    return coeffs, sol, ns


# ════════════════════════════════════════════════════════
# FEATURE 3: MULTI-SOLUTION HANDLER (NULLITY > 1)
# ════════════════════════════════════════════════════════

def solve_all(matrix):
    """
    Returns list of integer coefficient vectors — one per null space basis vector.
    Each is independently LCM+GCD reduced to minimal integers.
    """
    ns = matrix.nullspace()
    if not ns:
        return [], ns
    solutions = []
    for vec in ns:
        lcm_val = lcm([term.q for term in vec])
        coeffs  = [int(term * lcm_val) for term in vec]
        g = reduce(gcd, [abs(c) for c in coeffs])
        coeffs = [c // g for c in coeffs]
        solutions.append(coeffs)
    return solutions, ns


# ════════════════════════════════════════════════════════
# CORE HELPERS
# ════════════════════════════════════════════════════════

def build_matrix(reactants, products):
    compounds = reactants + products
    elements  = set()
    for c in compounds:
        elements.update(parse_compound(c).keys())
    elements = sorted(elements)
    matrix   = []
    for el in elements:
        row = []
        for c in reactants:
            row.append( parse_compound(c).get(el, 0))
        for c in products:
            row.append(-parse_compound(c).get(el, 0))
        matrix.append(row)
    return Matrix(matrix), elements


def check_elements(reactants, products):
    """Returns (balanced_bool, left_only_set, right_only_set)"""
    left  = set()
    right = set()
    for c in reactants:
        left.update(parse_compound(c).keys())
    for c in products:
        right.update(parse_compound(c).keys())
    return (left == right), (left - right), (right - left)


def format_eq(reactants, products, coeffs):
    r_len = len(reactants)
    left  = " + ".join(
        (f"{coeffs[i]}" if coeffs[i] != 1 else "") + reactants[i]
        for i in range(r_len)
    )
    right = " + ".join(
        (f"{coeffs[i + r_len]}" if coeffs[i + r_len] != 1 else "") + products[i]
        for i in range(len(products))
    )
    return left + "  →  " + right


def verify_balance(reactants, products, coeffs):
    """Numerically verify Ax = 0 — returns True if perfectly balanced."""
    all_els = set()
    for c in reactants + products:
        all_els.update(parse_compound(c).keys())
    for el in all_els:
        lhs = sum(coeffs[i] * parse_compound(r).get(el, 0) for i, r in enumerate(reactants))
        rhs = sum(coeffs[len(reactants)+i] * parse_compound(p).get(el, 0) for i, p in enumerate(products))
        if lhs != rhs:
            return False
    return True


# ════════════════════════════════════════════════════════
# SESSION STATE
# ════════════════════════════════════════════════════════

for key, default in [("run", False), ("chosen_sol", 0), ("last_eq", "")]:
    if key not in st.session_state:
        st.session_state[key] = default

# ════════════════════════════════════════════════════════
# UI — INPUT
# ════════════════════════════════════════════════════════

st.markdown(
    "Now supports **parentheses** · **minimal coefficients** · **multiple solutions** · **precise error reporting**"
)

col_input, col_btn = st.columns([5, 1])
with col_input:
    equation = st.text_input("equation", placeholder="e.g. Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O",
                             label_visibility="collapsed")
with col_btn:
    if st.button("Balance ⚗️", type="primary", use_container_width=True):
        st.session_state.run = True
        st.session_state.chosen_sol = 0
        st.session_state.last_eq = equation

st.caption("Quick examples:")
ex_cols = st.columns(5)
examples = [
    ("Ethane combustion",  "C2H6 + O2 -> CO2 + H2O"),
    ("Rust formation",     "Fe + O2 -> Fe2O3"),
    ("Calcium phosphate",  "Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O"),
    ("Aluminium sulfate",  "Al2(SO4)3 + NaOH -> Al(OH)3 + Na2SO4"),
    ("Photosynthesis",     "CO2 + H2O -> C6H12O6 + O2"),
]
for col, (label, eq) in zip(ex_cols, examples):
    if col.button(label, use_container_width=True):
        st.session_state.run = True
        st.session_state.chosen_sol = 0
        st.session_state.last_eq = eq
        equation = eq

if st.session_state.last_eq:
    equation = st.session_state.last_eq

# ════════════════════════════════════════════════════════
# MAIN ANALYSIS
# ════════════════════════════════════════════════════════

if st.session_state.run and equation:
    try:
        if "->" not in equation:
            st.error("Format:  Reactants -> Products")
            st.stop()

        left_str, right_str = equation.split("->", 1)
        reactants = [x.strip() for x in left_str.split("+") if x.strip()]
        products  = [x.strip() for x in right_str.split("+") if x.strip()]

        if not reactants or not products:
            st.error("Need at least one reactant and one product.")
            st.stop()

        # ── FEATURE 4: Precise element mismatch reporting ──
        ok, left_only, right_only = check_elements(reactants, products)
        if not ok:
            msg = "❌ Element mismatch — equation cannot be balanced.\n\n"
            if left_only:
                msg += f"**Only on LEFT side:** {', '.join(sorted(left_only))}\n\n"
            if right_only:
                msg += f"**Only on RIGHT side:** {', '.join(sorted(right_only))}"
            st.error(msg)
            st.stop()

        # ── Section 1: Parsed compounds (FEATURE 1 on display) ──
        st.subheader("🔷 Parsed Compounds")
        st.caption("Recursive parser correctly expands brackets like Ca(OH)₂ and Al₂(SO₄)₃.")
        p_cols = st.columns(min(len(reactants + products), 4))
        for idx, comp in enumerate(reactants + products):
            parsed = parse_compound(comp)
            role   = "Reactant" if idx < len(reactants) else "Product"
            bg     = "#e8f5e9" if role == "Reactant" else "#ede7f6"
            fg     = "#1b5e20" if role == "Reactant" else "#311b92"
            with p_cols[idx % len(p_cols)]:
                st.markdown(
                    f"<div style='border:1px solid #ccc;border-radius:8px;"
                    f"padding:10px;margin:4px 0;background:{bg}'>"
                    f"<b style='color:{fg}'>{comp}</b><br>"
                    f"<span style='font-size:11px;color:#555'>{role}</span><br>"
                    f"<code style='font-size:11px'>{dict(parsed)}</code></div>",
                    unsafe_allow_html=True
                )

        # ── Section 2: Matrix ──
        matrix, elements = build_matrix(reactants, products)
        n = len(reactants) + len(products)

        st.subheader("🔷 Element Conservation Matrix  (A)")
        st.caption("Rows = elements · Columns = compounds · Reactants → +ve · Products → −ve")
        df = pd.DataFrame(matrix.tolist(), index=elements, columns=reactants + products)
        st.dataframe(df, use_container_width=True)

        # ── Section 3: System of equations ──
        st.subheader("🔷 System of Linear Equations  (Ax = 0)")
        for i, el in enumerate(elements):
            terms = []
            for j, comp in enumerate(reactants + products):
                c = matrix[i, j]
                if c != 0:
                    terms.append(f"({c})·x{j+1}")
            st.write(f"**{el}:**   " + " + ".join(terms) + " = 0")

        # ── Section 4: RREF ──
        st.subheader("🔷 RREF  (Reduced Row Echelon Form)")
        rref_matrix, pivots = matrix.rref()
        st.write(rref_matrix)
        free_vars = [i for i in range(n) if i not in pivots]
        st.info(
            f"Pivot columns (determined): {list(pivots)}  ·  "
            f"Free variable columns: {free_vars}"
        )

        # ── Section 5: Rank & Nullity ──
        rank    = matrix.rank()
        nullity = n - rank
        st.subheader("🔷 Rank & Nullity")
        m1, m2, m3 = st.columns(3)
        m1.metric("Rank",      rank,    help="Independent element constraints")
        m2.metric("Nullity",   nullity, help="Degrees of freedom / independent solutions")
        m3.metric("Compounds", n)

        if nullity == 0:
            st.error("Nullity = 0 → no non-trivial solution. Cannot balance this equation.")
            st.stop()

        # ── Section 6: Null space ──
        st.subheader("🔷 Null Space")
        ns_vecs = matrix.nullspace()
        st.write(ns_vecs)
        st.caption(f"{nullity} basis vector(s). Each encodes one independent coefficient ratio.")

        # ── FEATURE 3: Multi-solution handler ──
        solutions, ns_raw = solve_all(matrix)
        if not solutions:
            st.error("No valid solution found.")
            st.stop()

        if nullity > 1:
            st.subheader(f"⚠️ Multiple Solutions  (Nullity = {nullity})")
            st.markdown(
                f"This equation has **{nullity} independent solution directions**. "
                "Select which basis vector to use:"
            )
            sol_labels = [
                f"Solution {i+1}: {format_eq(reactants, products, s)}"
                for i, s in enumerate(solutions)
            ]
            chosen = st.radio("Pick a solution:", range(len(solutions)),
                              format_func=lambda i: sol_labels[i])
            st.session_state.chosen_sol = chosen
            st.info(
                "Nullity > 1 often means the equation represents overlapping "
                "independent reactions. Each basis vector is one minimal sub-reaction."
            )
            coeffs = solutions[chosen]
        else:
            coeffs = solutions[0]

        # ── Section 7: Coefficient step-by-step (FEATURE 2) ──
        st.subheader("🔷 LCM Scale  →  GCD Reduce  →  Minimal Coefficients")
        raw_vec = ns_raw[st.session_state.chosen_sol if nullity > 1 else 0]
        denoms  = [term.q for term in raw_vec]
        lv      = lcm(denoms)
        pre_gcd = [int(term * lv) for term in raw_vec]
        g       = reduce(gcd, [abs(c) for c in pre_gcd])

        step_df = pd.DataFrame({
            "Compound":        reactants + products,
            "Null vector":     [str(t) for t in raw_vec],
            f"× LCM({lv})":   pre_gcd,
            f"÷ GCD({g})":    coeffs,
        })
        st.dataframe(step_df, use_container_width=True, hide_index=True)

        # ── Section 8: Balanced equation ──
        balanced = format_eq(reactants, products, coeffs)
        verified = verify_balance(reactants, products, coeffs)

        st.subheader("✅ Balanced Equation")
        st.success(balanced)
        if verified:
            st.caption("✔ Numerically verified: atom counts match on both sides (Ax = 0 confirmed).")
        else:
            st.warning("⚠ Verification failed — check input.")

        # ── Section 9: Mole Ratio Chart ──
        st.subheader("📊 Mole Ratio Visualisation")
        fig, ax = plt.subplots(figsize=(max(6, n * 1.3), 4))
        ax.set_facecolor("#0D1B2A")
        fig.patch.set_facecolor("#0D1B2A")
        colors = ["#00D9A3" if i < len(reactants) else "#7C3AED" for i in range(n)]
        bars = ax.bar(reactants + products, coeffs, color=colors, edgecolor="#1A2F4A", width=0.6)
        ax.set_ylabel("Stoichiometric coefficient", color="white")
        ax.tick_params(colors="white")
        for sp in ax.spines.values(): sp.set_color("#334D6E")
        for bar, val in zip(bars, coeffs):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                    str(val), ha="center", va="bottom", color="white", fontsize=10, fontweight="bold")
        ax.legend(handles=[
            mpatches.Patch(facecolor="#00D9A3", label="Reactants"),
            mpatches.Patch(facecolor="#7C3AED", label="Products"),
        ], facecolor="#1A2F4A", labelcolor="white", edgecolor="#334D6E")
        st.pyplot(fig)

        # ── Section 10: Limiting Reagent Simulator ──
        st.subheader("🧪 Limiting Reagent Simulator")
        st.caption("Enter moles of each reactant. App shows which runs out first and how much product forms.")

        user_inputs = []
        mol_cols = st.columns(len(reactants))
        for i, r in enumerate(reactants):
            with mol_cols[i]:
                val = st.number_input(f"Moles of {r}", min_value=0.001, value=1.0,
                                      step=0.1, key=f"mol_{i}_{r}")
                user_inputs.append(val)

        ratios  = [user_inputs[i] / coeffs[i] for i in range(len(reactants))]
        lim_idx = ratios.index(min(ratios))
        scale   = ratios[lim_idx]

        rc1, rc2 = st.columns(2)
        with rc1:
            st.error(f"**Limiting reagent:** `{reactants[lim_idx]}`")
            st.caption(f"Ratio vector (available ÷ coefficient): {[round(r, 4) for r in ratios]}")
            st.caption(f"Minimum ratio = {round(scale, 4)} → `{reactants[lim_idx]}` exhausted first")
        with rc2:
            st.markdown("**Products formed:**")
            for i, p in enumerate(products):
                st.write(f"→ `{p}`: **{scale * coeffs[len(reactants)+i]:.4f}** moles")
            st.markdown("**Excess reactants remaining:**")
            for i, r in enumerate(reactants):
                if i != lim_idx:
                    remaining = user_inputs[i] - scale * coeffs[i]
                    st.write(f"→ `{r}`: **{remaining:.4f}** moles unused")

    except Exception as e:
        st.error(f"Error: {e}")
        st.caption("Expected format: `Reactant1 + Reactant2 -> Product1 + Product2`")
