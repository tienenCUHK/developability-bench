#!/usr/bin/env python3
"""Static sanity checker for the DevMPNN-Bench LaTeX manuscript.

There is no TeX installation in this environment (no pdflatex/latexmk/pandoc),
so this script substitutes a compile with structural checks that catch the
mistakes a real compile would reject or warn about:

  1. environment balance      -- every \\begin{env} has a matching \\end{env}
  2. brace balance            -- { } balanced, never negative, escaped \\{ \\}
                                 and comments ignored
  3. math delimiter balance   -- $ ... $ and \\[ ... \\] and \\left/\\right
  4. cross-references         -- every \\ref/\\eqref/\\autoref/... target has a
                                 matching \\label; \\label names are unique
  5. figures                  -- every \\includegraphics path resolves on disk
                                 (honouring \\graphicspath, trying common
                                 extensions)
  6. inputs                   -- every \\input/\\include file exists
  7. bibliography             -- when refs.bib is present: every \\cite key is
                                 defined, expected keys are present, no
                                 duplicate/dead keys, and only the sanctioned
                                 citation keys are used
  8. encoding                 -- no non-ASCII bytes in .tex sources (keeps
                                 pdflatex + [utf8]{inputenc} safe)

Run from anywhere:   python3 check_latex.py [PAPER_DIR]
Exit code 0 = all checks pass (warnings allowed), 1 = at least one error.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Citation keys the task specification allows (see CITATION_MAP.md).
EXPECTED_CITE_KEYS = [
    "struc_transformer", "gvp", "proteinmpnn", "esm_if1", "pifold",
    "ligandmpnn", "soluble_mpnn", "hypermpnn", "membranempnn", "proteinguide",
    "discrete_guidance", "fk_rfdiffusion", "zetadial", "proteininvbench",
    "motifbench", "proteingym", "flip", "esmfold", "protenix",
]

REF_CMDS = ("ref", "eqref", "autoref", "cref", "Cref", "pageref", "nameref")
GRAPHICS_EXTS = (".png", ".pdf", ".jpg", ".jpeg", ".eps", "")

errors = []
warnings = []
notes = []


def err(msg):
    errors.append(msg)


def warn(msg):
    warnings.append(msg)


def rel(path):
    """Path relative to the paper directory, for readable messages."""
    try:
        return os.path.relpath(path, HERE)
    except ValueError:
        return path


def strip_comment(line):
    """Remove an unescaped % comment from one line (LaTeX rules)."""
    i = 0
    while i < len(line):
        if line[i] == "%":
            # Count immediately preceding backslashes: odd -> escaped percent.
            n = 0
            j = i - 1
            while j >= 0 and line[j] == "\\":
                n += 1
                j -= 1
            if n % 2 == 0:
                return line[:i]
        i += 1
    return line


class TexFile(object):
    """A comment-stripped .tex file with a char-offset -> line-number map."""

    def __init__(self, path):
        self.path = path
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            raw = fh.read()
        self.raw = raw
        stripped_lines = [strip_comment(ln) for ln in raw.splitlines()]
        self.text = "\n".join(stripped_lines)
        self.line_starts = [0]
        for ln in stripped_lines:
            self.line_starts.append(self.line_starts[-1] + len(ln) + 1)

    def line_of(self, offset):
        lo, hi = 0, len(self.line_starts) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.line_starts[mid] <= offset:
                lo = mid
            else:
                hi = mid - 1
        return lo + 1


def tex_files(main_path):
    """Return the main file plus every \\input-ed file, in document order."""
    out = []
    seen = set()

    def visit(path):
        path = os.path.abspath(path)
        if path in seen or not os.path.isfile(path):
            return
        seen.add(path)
        tf = TexFile(path)
        out.append(tf)
        for m in re.finditer(r"\\(?:input|include)\{([^}]*)\}", tf.text):
            target = m.group(1).strip()
            if not target:
                continue
            candidate = os.path.join(os.path.dirname(path), target)
            if not os.path.isfile(candidate) and not os.path.splitext(candidate)[1]:
                candidate += ".tex"
            if not os.path.isfile(candidate):
                err("%s:%d: \\input file not found: %s"
                    % (rel(path), tf.line_of(m.start()), target))
                continue
            visit(candidate)

    visit(main_path)
    return out


def check_document_wrapper(main_tf):
    if not re.search(r"\\documentclass", main_tf.text):
        err("main.tex: no \\documentclass found")
    if "\\begin{document}" not in main_tf.text:
        err("main.tex: no \\begin{document} found")
    if "\\end{document}" not in main_tf.text:
        err("main.tex: no \\end{document} found")
    if "\\begin{document}" in main_tf.text and "\\end{document}" in main_tf.text:
        if main_tf.text.index("\\begin{document}") > main_tf.text.index("\\end{document}"):
            err("main.tex: \\end{document} precedes \\begin{document}")
    for env in ("document",):
        nbeg = len(re.findall(r"\\begin\{%s\}" % env, main_tf.text))
        nend = len(re.findall(r"\\end\{%s\}" % env, main_tf.text))
        if nbeg != 1 or nend != 1:
            err("main.tex: expected exactly one %s environment, found %d begin / %d end"
                % (env, nbeg, nend))


def check_environment_balance(files):
    for tf in files:
        stack = []
        for m in re.finditer(r"\\(begin|end)\{([^}]*)\}", tf.text):
            kind, env = m.group(1), m.group(2).strip()
            line = tf.line_of(m.start())
            if kind == "begin":
                stack.append((env, line))
            else:
                if not stack:
                    err("%s:%d: \\end{%s} without matching \\begin"
                        % (rel(tf.path), line, env))
                elif stack[-1][0] != env:
                    open_env, open_line = stack[-1]
                    err("%s:%d: \\end{%s} closes \\begin{%s} opened at line %d"
                        % (rel(tf.path), line, env, open_env, open_line))
                    stack.pop()
                else:
                    stack.pop()
        for env, line in stack:
            err("%s:%d: \\begin{%s} is never closed"
                % (rel(tf.path), line, env))


def check_braces(files):
    for tf in files:
        depth = 0
        i = 0
        text = tf.text
        while i < len(text):
            c = text[i]
            if c == "\\":
                i += 2
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth < 0:
                    err("%s:%d: unbalanced '}' (no matching '{')"
                        % (rel(tf.path), tf.line_of(i)))
                    depth = 0
            i += 1
        if depth != 0:
            err("%s: %d unclosed '{' at end of file" % (rel(tf.path), depth))


def check_math(files):
    for tf in files:
        # Inline math: an even number of unescaped $ in each file.
        dollars = len(re.findall(r"(?<!\\)\$", tf.text))
        if dollars % 2:
            err("%s: odd number of unescaped '$' (%d) -- inline math is unbalanced"
                % (rel(tf.path), dollars))
        n_open = len(re.findall(r"(?<!\\)\\\[", tf.text))
        n_close = len(re.findall(r"(?<!\\)\\\]", tf.text))
        if n_open != n_close:
            err("%s: \\[ ... \\] unbalanced (%d open, %d close)"
                % (rel(tf.path), n_open, n_close))
        n_left = len(re.findall(r"\\left\b", tf.text))
        n_right = len(re.findall(r"\\right\b", tf.text))
        if n_left != n_right:
            err("%s: \\left/\\right unbalanced (%d \\left, %d \\right)"
                % (rel(tf.path), n_left, n_right))


def check_refs(files):
    labels = {}
    for tf in files:
        for m in re.finditer(r"\\label\{([^}]*)\}", tf.text):
            name = m.group(1).strip()
            line = tf.line_of(m.start())
            if name in labels:
                err("%s:%d: duplicate \\label{%s} (first defined at %s:%d)"
                    % (rel(tf.path), line, name, rel(labels[name][0]),
                       labels[name][1]))
            else:
                labels[name] = (tf.path, line)
    for tf in files:
        pat = r"\\(%s)\*?(?:\[[^\]]*\])?\{([^}]*)\}" % "|".join(REF_CMDS)
        for m in re.finditer(pat, tf.text):
            cmd, targets = m.group(1), m.group(2)
            for target in [t.strip() for t in targets.split(",") if t.strip()]:
                if target not in labels:
                    err("%s:%d: \\%s{%s} has no matching \\label"
                        % (rel(tf.path), tf.line_of(m.start()), cmd, target))


def parse_graphicspath(main_tf):
    paths = []
    for m in re.finditer(r"\\graphicspath\{((?:\{[^}]*\})+)\}", main_tf.text):
        paths.extend(re.findall(r"\{([^}]*)\}", m.group(1)))
    return paths


def check_graphics(files, main_tf):
    base = os.path.dirname(os.path.abspath(main_tf.path))
    graphicspaths = parse_graphicspath(main_tf)
    found = []
    pat = r"\\includegraphics\*?(?:\[[^\]]*\])?\{([^}]*)\}"
    for tf in files:
        for m in re.finditer(pat, tf.text):
            raw_path = m.group(1).strip()
            line = tf.line_of(m.start())
            roots = [base] + [os.path.normpath(os.path.join(base, gp))
                              for gp in graphicspaths]
            resolved = None
            for root in roots:
                cand = os.path.normpath(os.path.join(root, raw_path))
                for ext in GRAPHICS_EXTS:
                    if ext and os.path.splitext(cand)[1]:
                        trial = cand
                    elif ext:
                        trial = cand + ext
                    else:
                        trial = cand
                    if os.path.isfile(trial):
                        resolved = trial
                        break
                if resolved:
                    break
            if resolved is None:
                err("%s:%d: \\includegraphics{%s} does not resolve to an existing file"
                    % (rel(tf.path), line, raw_path))
            else:
                found.append(resolved)
    # Every figure shipped in paper/figs should be used at least once.
    fig_dir = os.path.join(base, "figs")
    if os.path.isdir(fig_dir):
        for name in sorted(os.listdir(fig_dir)):
            full = os.path.join(fig_dir, name)
            if os.path.isfile(full) and os.path.abspath(full) not in found:
                warn("figure not referenced by any \\includegraphics: %s" % rel(full))


def parse_bib_keys(bib_path):
    with open(bib_path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    keys = []
    for m in re.finditer(r"@\s*([A-Za-z]+)\s*[{(]\s*([^,\s}]+)\s*,", text):
        entry_type, key = m.group(1).lower(), m.group(2)
        if entry_type in ("string", "comment", "preamble"):
            continue
        keys.append(key)
    return text, keys


def check_bibliography(files, main_tf):
    base = os.path.dirname(os.path.abspath(main_tf.path))
    bib_name = "refs"
    m = re.search(r"\\bibliography\{([^}]*)\}", main_tf.text)
    if m:
        bib_name = m.group(1).split(",")[0].strip()
    bib_path = os.path.join(base, bib_name + ".bib")
    bib_style = re.search(r"\\bibliographystyle\{([^}]*)\}", main_tf.text)
    if not bib_style:
        warn("main.tex: no \\bibliographystyle found")
    if not os.path.isfile(bib_path):
        warn("bibliography not found yet: %s (run again once refs.bib exists)"
             % rel(bib_path))
        return
    text, bib_keys = parse_bib_keys(bib_path)
    seen = set()
    for key in bib_keys:
        if key in seen:
            err("%s: duplicate bibliography key '%s'" % (rel(bib_path), key))
        seen.add(key)
    missing_expected = [k for k in EXPECTED_CITE_KEYS if k not in seen]
    if missing_expected:
        warn("%s: expected keys not defined yet: %s"
             % (rel(bib_path), ", ".join(missing_expected)))
    unbraced = 0
    depth = 0
    i = 0
    while i < len(text):
        c = text[i]
        if c == "\\":
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            unbraced = min(unbraced, depth)
        i += 1
    if depth != 0 or unbraced < 0:
        err("%s: unbalanced braces in bibliography" % rel(bib_path))
    for tf in files:
        for m in re.finditer(r"\\cite[a-zA-Z]*\*?(?:\[[^\]]*\])*\{([^}]*)\}",
                             tf.text):
            line = tf.line_of(m.start())
            for key in [k.strip() for k in m.group(1).split(",") if k.strip()]:
                if key not in seen:
                    err("%s:%d: \\cite{%s} is not defined in %s"
                        % (rel(tf.path), line, key, os.path.basename(bib_path)))
                if key not in EXPECTED_CITE_KEYS:
                    err("%s:%d: \\cite{%s} is not one of the sanctioned citation keys"
                        % (rel(tf.path), line, key))
    cited = set()
    for tf in files:
        for m in re.finditer(r"\\cite[a-zA-Z]*\*?(?:\[[^\]]*\])*\{([^}]*)\}",
                             tf.text):
            cited.update(k.strip() for k in m.group(1).split(",") if k.strip())
    unused = [k for k in EXPECTED_CITE_KEYS if k in seen and k not in cited]
    if unused:
        warn("sanctioned citation keys defined but never cited: %s"
             % ", ".join(unused))
    notes.append("bibliography: %d entries, %d distinct citation keys used"
                 % (len(seen), len(cited)))


def check_ascii(files):
    for tf in files:
        with open(tf.path, "rb") as fh:
            data = fh.read()
        bad_lines = []
        for lineno, raw_line in enumerate(data.split(b"\n"), 1):
            if any(b > 127 for b in raw_line):
                bad_lines.append(lineno)
        if bad_lines:
            shown = ", ".join(str(n) for n in bad_lines[:8])
            extra = "" if len(bad_lines) <= 8 else " ..."
            err("%s: non-ASCII bytes on line(s) %s%s (pdflatex + inputenc "
                "utf8 is safest with pure-ASCII sources)" % (rel(tf.path), shown, extra))


def main():
    paper_dir = sys.argv[1] if len(sys.argv) > 1 else HERE
    paper_dir = os.path.abspath(paper_dir)
    main_path = os.path.join(paper_dir, "main.tex")
    if not os.path.isfile(main_path):
        sys.stderr.write("ERROR: %s not found\n" % main_path)
        return 1

    main_tf = TexFile(main_path)
    files = tex_files(main_path)

    check_document_wrapper(main_tf)
    check_environment_balance(files)
    check_braces(files)
    check_math(files)
    check_refs(files)
    check_graphics(files, main_tf)
    check_bibliography(files, main_tf)
    check_ascii(files)

    print("=" * 72)
    print("check_latex.py -- static checks for %s" % rel(main_path))
    print("=" * 72)
    print("files checked (%d):" % len(files))
    for tf in files:
        print("  %s" % rel(tf.path))
    print("-" * 72)
    if notes:
        for n in notes:
            print("NOTE   %s" % n)
    if warnings:
        for w in warnings:
            print("WARN   %s" % w)
    if errors:
        for e in errors:
            print("ERROR  %s" % e)
        print("-" * 72)
        print("RESULT: FAIL (%d error(s), %d warning(s))"
              % (len(errors), len(warnings)))
        return 1
    print("RESULT: PASS -- environment/brace/math balance OK, all \\ref "
          "targets labeled,")
    print("        all graphics resolved, all \\input files present, "
          "citation keys consistent.")
    if warnings:
        print("        (%d warning(s) above are advisory)" % len(warnings))
    return 0


if __name__ == "__main__":
    sys.exit(main())
