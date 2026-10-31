const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, ImageRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  LineNumberRestartFormat, convertInchesToTwip, Footer, PageNumber,
} = require("docx");

const SRC = process.argv[2] || "manuscript.md";
const OUT = process.argv[3] || "manuscript.docx";
const FONT = "Times New Roman";
const PAGE_W = 12240, PAGE_H = 15840, MARGIN = 1440;       // US Letter, 1 inch
const USABLE_IN = (PAGE_W - 2 * MARGIN) / 1440;             // 6.5 inches
const dims = JSON.parse(fs.readFileSync("figs.json", "utf8"));

/* ---------- inline markdown -> TextRuns ---------- */
function runs(text, base = {}) {
  const out = [];
  // ^a^ superscript, **bold**, *italic*
  text = text.replace(/\\([*^_])/g, "$1").replace(/`([^`]+)`/g, "$1");
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]+\^)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) out.push(new TextRun({ ...base, text: text.slice(last, m.index), font: FONT }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new TextRun({ ...base, text: t.slice(2, -2), bold: true, font: FONT }));
    else if (t.startsWith("^")) out.push(new TextRun({ ...base, text: t.slice(1, -1), superScript: true, font: FONT }));
    else out.push(new TextRun({ ...base, text: t.slice(1, -1), italics: true, font: FONT }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({ ...base, text: text.slice(last), font: FONT }));
  return out.length ? out : [new TextRun({ ...base, text: "", font: FONT })];
}

const para = (text, opts = {}) => new Paragraph({
  children: runs(text),
  spacing: { line: opts.line ?? 480, after: opts.after ?? 120 },   // 480 = double
  alignment: opts.alignment,
  indent: opts.indent,
  heading: opts.heading,
  keepNext: opts.keepNext,
  pageBreakBefore: opts.pageBreakBefore,
});

const heading = (text, level) => new Paragraph({
  children: runs(text, { bold: true, size: level === HeadingLevel.HEADING_1 ? 26 : 24 }),
  heading: level,
  spacing: { before: 300, after: 160, line: 360 },
  keepNext: true,
});

/* ---------- build ---------- */
const md = fs.readFileSync(SRC, "utf8").split(/\r?\n/);
const children = [];
let i = 0;

while (i < md.length) {
  const line = md[i];

  if (/^---\s*$/.test(line)) { i++; continue; }
  if (/^\s*$/.test(line)) { i++; continue; }

  // image
  let im = line.match(/^!\[.*?\]\((.+?)\)\s*$/);
  if (im) {
    const file = im[1];
    const base = path.basename(file);
    const d = dims[base];
    const wIn = Math.min(USABLE_IN, 6.5);
    const hIn = d ? wIn * (d.h / d.w) : 3;
    children.push(new Paragraph({
      children: [new ImageRun({
        type: "png",
        data: fs.readFileSync(file),
        transformation: { width: Math.round(wIn * 96), height: Math.round(hIn * 96) },
      })],
      alignment: AlignmentType.CENTER,
      spacing: { before: 200, after: 80 },
      keepNext: true,
    }));
    i++; continue;
  }

  // table
  if (/^\|/.test(line) && i + 1 < md.length && /^\|[\s:|-]+\|/.test(md[i + 1])) {
    const rows = [];
    while (i < md.length && /^\|/.test(md[i])) {
      if (!/^\|[\s:|-]+\|$/.test(md[i].trim()))
        rows.push(md[i].trim().replace(/^\||\|$/g, "").split("|").map(c => c.trim()));
      i++;
    }
    const ncol = Math.max(...rows.map(r => r.length));
    const tw = PAGE_W - 2 * MARGIN;
    const colw = Array(ncol).fill(Math.floor(tw / ncol));
    colw[0] = tw - colw.slice(1).reduce((a, b) => a + b, 0);
    children.push(new Table({
      columnWidths: colw,
      width: { size: tw, type: WidthType.DXA },
      rows: rows.map((r, ri) => new TableRow({
        tableHeader: ri === 0,
        cantSplit: true,
        children: Array.from({ length: ncol }, (_, ci) => new TableCell({
          width: { size: colw[ci], type: WidthType.DXA },
          shading: ri === 0 ? { type: ShadingType.CLEAR, fill: "F2F2F2" } : undefined,
          margins: { top: 60, bottom: 60, left: 100, right: 100 },
          children: [new Paragraph({
            children: runs(r[ci] ?? "", { bold: ri === 0, size: 20 }),
            spacing: { line: 240, after: 0 },
          })],
        })),
      })),
    }));
    children.push(new Paragraph({ text: "", spacing: { after: 160 } }));
    continue;
  }

  // headings
  if (/^# /.test(line)) {
    children.push(new Paragraph({
      children: runs(line.slice(2), { bold: true, size: 30 }),
      spacing: { after: 240, line: 480 },
      alignment: AlignmentType.LEFT,
    }));
    i++; continue;
  }
  if (/^## /.test(line)) { children.push(heading(line.slice(3), HeadingLevel.HEADING_1)); i++; continue; }
  if (/^### /.test(line)) { children.push(heading(line.slice(4), HeadingLevel.HEADING_2)); i++; continue; }

  // bullet
  if (/^- /.test(line)) {
    children.push(new Paragraph({
      children: runs(line.slice(2)),
      bullet: { level: 0 },
      spacing: { line: 360, after: 60 },
    }));
    i++; continue;
  }

  // figure / table caption
  if (/^\*\*(Fig\.|Table)/.test(line)) {
    children.push(new Paragraph({
      children: runs(line, { size: 20 }),
      spacing: { line: 480, after: 200 },
      keepNext: /^\*\*Table/.test(line),
    }));
    i++; continue;
  }

  // reference entry: hanging indent
  const inRefs = children.some(c => c instanceof Paragraph) &&
    md.slice(0, i).some(l => /^## References/.test(l));
  if (inRefs && /^[a-zA-ZÇÁÉÍÓÚÀÈÌÒÙŠŽČĆĐ]/.test(line) && /\((?:19|20)\d{2}[a-z]?\)/.test(line)) {
    children.push(new Paragraph({
      children: runs(line, { size: 20 }),
      spacing: { line: 480, after: 120 },
      indent: { left: convertInchesToTwip(0.4), hanging: convertInchesToTwip(0.4) },
    }));
    i++; continue;
  }

  children.push(para(line));
  i++;
}

const doc = new Document({
  styles: {
    default: {
      document: { run: { font: FONT, size: 24, color: "000000" } },
      heading1: { run: { font: FONT, size: 26, bold: true, color: "000000" } },
      heading2: { run: { font: FONT, size: 24, bold: true, color: "000000" } },
    },
  },
  sections: [{
    properties: {
      page: { size: { width: PAGE_W, height: PAGE_H }, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } },
      lineNumbers: { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS, distance: 360 },
    },
    footers: {
      default: new Footer({
        children: [new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })],
        })],
      }),
    },
    children,
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync(OUT, b);
  console.log(`wrote ${OUT} (${(b.length / 1024).toFixed(0)} KB, ${children.length} blocks)`);
});
