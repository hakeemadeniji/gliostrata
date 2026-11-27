"""Render the OSF supplement documents (docs/osf_supplements/*.md) to PDF.

A deliberately small Markdown subset: headings, paragraphs, bullet and numbered lists,
pipe tables, **bold**, *italic* and `code`. It is enough for the three supplements and has
no dependency beyond PyMuPDF.

    python scripts/render_osf_pdfs.py --out ../gliostrata_osf_bundle
"""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

import fitz  # PyMuPDF

DOCS = {
    "docs/osf_supplements/measured_variables.md": "GlioStrata_OSF_measured_variables.pdf",
    "docs/osf_supplements/indices.md": "GlioStrata_OSF_indices.pdf",
    "docs/osf_supplements/statistical_models.md": "GlioStrata_OSF_statistical_models.pdf",
}
CSS = """
body { font-family: sans-serif; font-size: 9pt; line-height: 1.35; }
h1 { font-size: 15pt; margin-bottom: 6pt; }
h2 { font-size: 12pt; margin-top: 10pt; margin-bottom: 4pt; }
h3 { font-size: 10pt; margin-top: 8pt; }
p { margin-top: 3pt; margin-bottom: 3pt; }
table { border-collapse: collapse; margin-top: 4pt; margin-bottom: 6pt; }
th, td { border: 0.5pt solid #888; padding: 2pt 3pt; font-size: 7.5pt; vertical-align: top; }
th { background-color: #eeeeee; }
code { font-family: monospace; font-size: 8pt; }
li { margin-bottom: 2pt; }
"""


def _inline(text: str) -> str:
    out = html.escape(text, quote=False)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    out = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", out)
    out = re.sub(r"(?<![*\w])\*([^*]+)\*(?!\w)", r"<i>\1</i>", out)
    return out


def md_to_html(md: str) -> str:
    lines = md.splitlines()
    parts: list[str] = []
    para: list[str] = []
    i = 0

    def flush():
        if para:
            parts.append("<p>" + _inline(" ".join(para)) + "</p>")
            para.clear()

    while i < len(lines):
        line = lines[i]
        s = line.strip()
        if not s:
            flush()
            i += 1
            continue
        m = re.match(r"^(#{1,3})\s+(.*)$", s)
        if m:
            flush()
            n = len(m.group(1))
            parts.append(f"<h{n}>{_inline(m.group(2))}</h{n}>")
            i += 1
            continue
        if s.startswith("|"):
            flush()
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(lines[i].strip())
                i += 1
            cells = [[c.strip() for c in r.strip("|").split("|")] for r in rows]
            body = [r for r in cells if not all(re.fullmatch(r":?-{3,}:?", c) for c in r)]
            out = ["<table>"]
            for k, r in enumerate(body):
                tag = "th" if k == 0 else "td"
                out.append("<tr>" + "".join(f"<{tag}>{_inline(c)}</{tag}>" for c in r) + "</tr>")
            out.append("</table>")
            parts.append("".join(out))
            continue
        if re.match(r"^(- |\d+\. )", s):
            flush()
            # (indent, ordered, text) per item; a non-marker line continues the last item.
            items: list[list] = []
            while i < len(lines) and lines[i].strip():
                raw = lines[i]
                cur = raw.strip()
                if re.match(r"^(- |\d+\. )", cur):
                    indent = len(raw) - len(raw.lstrip(" "))
                    items.append([indent, bool(re.match(r"^\d+\. ", cur)),
                                  re.sub(r"^(- |\d+\. )", "", cur)])
                else:
                    items[-1][2] += " " + cur
                i += 1
            out: list[str] = []
            stack: list[tuple[int, str]] = []
            for indent, ordered, text in items:
                while stack and indent < stack[-1][0]:
                    out.append(f"</li></{stack.pop()[1]}>")
                if not stack or indent > stack[-1][0]:
                    tag = "ol" if ordered else "ul"
                    stack.append((indent, tag))
                    out.append(f"<{tag}>")
                else:
                    out.append("</li>")
                out.append(f"<li>{_inline(text)}")
            while stack:
                out.append(f"</li></{stack.pop()[1]}>")
            parts.append("".join(out))
            continue
        para.append(s)
        i += 1
    flush()
    return "<html><body>" + "".join(parts) + "</body></html>"


def render(md_path: Path, pdf_path: Path) -> int:
    story = fitz.Story(html=md_to_html(md_path.read_text(encoding="utf-8")), user_css=CSS)
    writer = fitz.DocumentWriter(str(pdf_path))
    mediabox = fitz.paper_rect("letter")
    where = mediabox + (54, 54, -54, -54)
    pages = 0
    more = True
    while more:
        device = writer.begin_page(mediabox)
        more, _ = story.place(where)
        story.draw(device)
        writer.end_page()
        pages += 1
    writer.close()
    return pages


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    for src, name in DOCS.items():
        n = render(Path(src), args.out / name)
        print(f"{name}: {n} page(s)")


if __name__ == "__main__":
    main()
