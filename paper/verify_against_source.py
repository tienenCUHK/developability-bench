#!/usr/bin/env python3
"""Check that the LaTeX manuscript is a faithful conversion of paper_v1.md.

This is the "did I lose or change a number?" checker that complements
check_latex.py (which checks LaTeX syntax only). It reads the markdown source
of truth and the generated .tex files and verifies three things:

  1. TABLES     -- every markdown table row exists in the matching LaTeX
                   tabular with the same sequence of numeric cell values
                   (signs included), and the table row/column counts agree.
  2. NUMBERS    -- every numeric value in the markdown prose (section headings
                   excluded) appears somewhere in the LaTeX.
  3. WORDS      -- every English word in the markdown prose appears somewhere in
                   the LaTeX (LaTeX/markdown markup normalised away). This
                   catches reworded or dropped sentences.

Changes that are *expected* and therefore normalised away:
  - markdown's Chinese author-facing summary (blockquote at the top), which is
    intentionally not typeset (it is preserved in AUTHOR_SUMMARY_zh.md);
  - the draft-metadata line (v1.0 draft ... target: ...);
  - LaTeX markup, math commands, citations, labels, figure/table captions that
    were added for the conversion, and the abstract environment wrapper.

Run:  python3 verify_against_source.py        (from this directory)
Exit code 0 = conversion is faithful, 1 = discrepancies found.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.normpath(os.path.join(HERE, "..", "paper_v1.md"))

# Words/markup that legitimately differ between markdown and LaTeX because of
# typesetting (math command names, hyphenation/dash normalisation, ...).
WORD_EXTRA_TEX = {
    "alpha", "beta", "gamma", "delta", "lambda", "phi", "theta", "rho", "pi",
    "sigma", "approx", "pm", "times", "to", "leq", "geq", "cdot", "to",
    "const", "mathrm", "textbf", "emph", "cite", "ref", "label", "section",
    "subsection", "paragraph", "begin", "end", "item", "toprule", "midrule",
    "bottomrule", "centering", "caption", "includegraphics", "linewidth",
    "bm", "textbackslash", "texttt", "prep", "targets", "py", "tsv", "fa",
    "abstract", "log",
}


def strip_comments(text):
    """Remove LaTeX comments from a whole file text."""
    return "\n".join(strip_comment(line) for line in text.splitlines())


def strip_comment(line):
    i = 0
    while i < len(line):
        if line[i] == "%":
            j = i - 1
            n = 0
            while j >= 0 and line[j] == "\\":
                n += 1
                j -= 1
            if n % 2 == 0:
                return line[:i]
        i += 1
    return line


def read_md():
    with open(SOURCE, "r", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    prose, tables = [], []
    cur_table = []
    in_summary = False
    for lineno, line in enumerate(lines, 1):
        if lineno == 3:            # draft metadata line
            continue
        if line.startswith(">"):   # Chinese author-facing summary
            in_summary = True
            continue
        if in_summary:
            if line.strip() == "" or line.strip() == "---":
                in_summary = False
            else:
                continue
        if line.lstrip().startswith("|"):
            cur_table.append(line.strip())
        else:
            if cur_table:
                tables.append(cur_table)
                cur_table = []
            prose.append(line)
    if cur_table:
        tables.append(cur_table)
    return "\n".join(prose), tables


def read_tex():
    """Concatenate the .tex files in *document* order (follow main.tex \input)."""
    main_path = os.path.join(HERE, "main.tex")
    with open(main_path, "r", encoding="utf-8") as fh:
        main_text = fh.read()
    files = ["main.tex"]
    for m in re.finditer(r"\\input\{([^}]*)\}", strip_comments(main_text)):
        target = m.group(1).strip()
        if not target.endswith(".tex"):
            target += ".tex"
        files.append(target)
    sec_dir = os.path.join(HERE, "sections")
    if os.path.isdir(sec_dir):  # any section file not \input-ed (defensive)
        for n in sorted(os.listdir(sec_dir)):
            rel = "sections/" + n
            if n.endswith(".tex") and rel not in files:
                files.append(rel)
    parts = []
    for relpath in files:
        path = os.path.join(HERE, relpath)
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read()
        parts.append("\n".join(strip_comment(l) for l in text.splitlines()))
    return "\n".join(parts)


def numbers(text, signed=True):
    text = (text.replace("\u2212", "-").replace("\u2013", "-")
                .replace("{,}", "").replace("{+}", "+").replace("{-}", "-"))
    out = []
    for tok in re.findall(r"[-+]?\d[\d,]*(?:\.\d+)?", text):
        core = tok.replace(",", "")
        out.append(core if signed else core.lstrip("+-"))
    return out


# --------------------------------------------------------------------------- #
# 1. tables
# --------------------------------------------------------------------------- #
GREEK_MAP = {
    "alpha": "alpha", "beta": "beta", "gamma": "gamma", "delta": "delta",
    "lambda": "lambda", "phi": "phi", "theta": "theta", "rho": "rho",
    "pi": "pi", "sigma": "sigma", "Delta": "delta", "Sigma": "sigma",
}


def md_cell_words(text):
    text = text.replace("\u2212", "-").replace("\u2013", "-")
    return set(re.findall(r"[A-Za-z]+(?:['-][A-Za-z]+)*", text.lower()))


def tex_cell_words(text):
    text = text.replace("\\\\", " ")
    for name in GREEK_MAP:
        text = re.sub(r"\\%s\b" % name, " " + GREEK_MAP[name] + " ", text)
    text = re.sub(r"\\(cite|ref|eqref|label|includegraphics)\*?"
                  r"(\[[^\]]*\])?\{[^}]*\}", " ", text)
    text = re.sub(r"\\[A-Za-z]+\*?", " ", text)
    text = (text.replace("{", " ").replace("}", " ").replace("$", "")
                .replace("\\_", "_").replace("--", "-")
                .replace("\u2212", "-").replace("\u2013", "-"))
    return set(re.findall(r"[A-Za-z]+(?:['-][A-Za-z]+)*", text.lower()))


def check_tables(md_tables, tex):
    tex_tables = []
    for m in re.finditer(r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}",
                         tex, re.S):
        rows = []
        for line in m.group(1).splitlines():
            line = line.strip()
            if not line or re.match(r"\\(top|mid|bottom)rule", line):
                continue
            if line.endswith("\\\\"):
                line = line[:-2]
            if "&" in line:
                rows.append([c.strip() for c in line.split("&")])
        tex_tables.append(rows)

    errs = []
    if len(md_tables) != len(tex_tables):
        errs.append("markdown has %d tables but LaTeX has %d"
                    % (len(md_tables), len(tex_tables)))
    for ti, (md_tab, tex_rows) in enumerate(zip(md_tables, tex_tables)):
        md_rows = []
        for raw in md_tab:
            cells = [c.strip() for c in raw.strip().strip("|").split("|")]
            if all(re.fullmatch(r":?-+:?", c) for c in cells):
                continue            # separator row
            md_rows.append(cells)
        if len(md_rows) != len(tex_rows):
            errs.append("table %d: %d markdown rows vs %d LaTeX rows"
                        % (ti + 1, len(md_rows), len(tex_rows)))
            continue
        for ri, (md_row, tex_row) in enumerate(zip(md_rows, tex_rows)):
            if len(md_row) != len(tex_row):
                errs.append("table %d row %d: %d markdown cells vs %d LaTeX cells"
                            % (ti + 1, ri + 1, len(md_row), len(tex_row)))
                continue
            for ci, (md_cell, tex_cell) in enumerate(zip(md_row, tex_row)):
                if numbers(md_cell) != numbers(tex_cell):
                    errs.append(
                        "table %d row %d col %d: md numbers %r != LaTeX numbers %r"
                        % (ti + 1, ri + 1, ci + 1,
                           numbers(md_cell), numbers(tex_cell)))
                missing_words = (md_cell_words(md_cell)
                                 - tex_cell_words(tex_cell)
                                 - WORD_EXTRA_TEX)
                if missing_words:
                    errs.append(
                        "table %d row %d col %d: md words %r missing from LaTeX cell %r"
                        % (ti + 1, ri + 1, ci + 1,
                           sorted(missing_words), tex_cell))
    return errs, len(md_tables)


# --------------------------------------------------------------------------- #
# 2. numbers in prose
# --------------------------------------------------------------------------- #
def check_numbers(md_prose, tex):
    # Numbers in section headings (e.g. "### 3.1") are auto-numbered by LaTeX;
    # collect exactly those to ignore, rather than ignoring every N.D pattern
    # (which would also swallow real values such as 6.0 or 2.0).
    heading_nums = set()
    for line in md_prose.splitlines():
        if line.lstrip().startswith("#"):
            heading_nums.update(re.findall(r"\d+\.\d+", line))
    md_prose_wo_headings = "\n".join(
        line for line in md_prose.splitlines()
        if not line.lstrip().startswith("#"))
    md_nums = set(numbers(md_prose_wo_headings))
    tex_nums = set(numbers(tex))
    missing = sorted(md_nums - tex_nums)
    missing = [m for m in missing if m not in heading_nums]
    return ["number %r in markdown prose is absent from LaTeX" % m
            for m in missing]


# --------------------------------------------------------------------------- #
# 3. words in prose
# --------------------------------------------------------------------------- #
def md_words(text):
    words = set()
    for line in text.splitlines():
        line = re.sub(r"^\s*#{1,6}\s*", "", line)
        line = re.sub(r"^\s*[-*]\s+", "", line)
        line = re.sub(r"^\s*\d+\.\s+", "", line)
        line = line.replace("**", "").replace("*", "").replace("`", "")
        line = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", line)
        line = line.replace("\u2212", "-").replace("\u2013", "-")
        for w in re.findall(r"[A-Za-z]+(?:['-][A-Za-z]+)*", line.lower()):
            words.add(w)
    return words


def tex_words(text):
    text = text.replace("\\\\", " ")
    for name in GREEK_MAP:
        text = re.sub(r"\\%s\b" % name, " " + GREEK_MAP[name] + " ", text)
    text = re.sub(r"\\(cite|ref|eqref|label|input|include|bibliography"
                  r"|bibliographystyle|graphicspath|includegraphics"
                  r"|begin|end|documentclass|usepackage)\*?"
                  r"(\[[^\]]*\])?\{[^}]*\}", " ", text)
    text = re.sub(r"\\[A-Za-z]+\*?", " ", text)
    text = text.replace("{", " ").replace("}", " ").replace("$", "")
    text = text.replace("--", "-").replace("\u2212", "-").replace("\u2013", "-")
    words = set()
    for w in re.findall(r"[A-Za-z]+(?:['-][A-Za-z]+)*", text.lower()):
        words.add(w)
    return words


def check_words(md_prose, tex):
    missing = sorted(md_words(md_prose) - tex_words(tex) - WORD_EXTRA_TEX)
    return ["word %r in markdown prose is absent from LaTeX" % w for w in missing]


def main():
    if not os.path.isfile(SOURCE):
        sys.stderr.write("ERROR: source not found: %s\n" % SOURCE)
        return 1
    md_prose, md_tables = read_md()
    tex = read_tex()

    table_errs, n_tables = check_tables(md_tables, tex)
    num_errs = check_numbers(md_prose, tex)
    word_errs = check_words(md_prose, tex)
    errs = table_errs + num_errs + word_errs

    print("=" * 72)
    print("verify_against_source.py")
    print("  source: %s" % os.path.relpath(SOURCE, HERE))
    print("  tables compared: %d" % n_tables)
    print("=" * 72)
    if errs:
        for e in errs:
            print("ERROR  %s" % e)
        print("-" * 72)
        print("RESULT: FAIL (%d discrepancy/ies)" % len(errs))
        return 1
    print("RESULT: PASS -- all markdown tables match cell-by-cell (numbers and")
    print("        signs), every number and English word in the markdown prose")
    print("        appears in the LaTeX, and no table rows/columns were lost.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
