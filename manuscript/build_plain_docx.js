/* A plain document: no line numbers, no page-number field, 1.15 spacing.
   Markdown understood: "# " title, "**bold**", "*italic*", "> " block quote,
   and a line of "key: value" metadata kept tight. */
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, AlignmentType, convertInchesToTwip,
} = require("docx");

const SRC = process.argv[2] || "response.md";
const OUT = process.argv[3] || "response.docx";
const FONT = "Times New Roman";
const LINE = 276;                       // 1.15 spacing

function runs(text, base = {}) {
  const out = [];
  text = text.replace(/\\([*^_])/g, "$1").replace(/`([^`]+)`/g, "$1");
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last)
      out.push(new TextRun({ ...base, text: text.slice(last, m.index), font: FONT }));
    const t = m[0];
    if (t.startsWith("**"))
      out.push(new TextRun({ ...base, text: t.slice(2, -2), bold: true, font: FONT }));
    else
      out.push(new TextRun({ ...base, text: t.slice(1, -1), italics: true, font: FONT }));
    last = m.index + t.length;
  }
  if (last < text.length)
    out.push(new TextRun({ ...base, text: text.slice(last), font: FONT }));
  return out.length ? out : [new TextRun({ ...base, text, font: FONT })];
}

const children = [];
for (const raw of fs.readFileSync(SRC, "utf8").split(/\r?\n/)) {
  const line = raw.trim();
  if (!line) continue;

  if (/^# /.test(line)) {
    children.push(new Paragraph({
      children: runs(line.slice(2), { bold: true, size: 26 }),
      spacing: { after: 240, line: LINE },
    }));
    continue;
  }
  if (/^> /.test(line)) {                       // quoted passage
    children.push(new Paragraph({
      children: runs(line.slice(2), { size: 22 }),
      indent: { left: convertInchesToTwip(0.4), right: convertInchesToTwip(0.4) },
      spacing: { before: 120, after: 200, line: LINE },
      alignment: AlignmentType.LEFT,
    }));
    continue;
  }
  if (/^[A-Za-z ]+:\s/.test(line) && line.length < 160) {   // metadata line
    children.push(new Paragraph({
      children: runs(line, { size: 22 }),
      spacing: { after: 40, line: LINE },
    }));
    continue;
  }
  children.push(new Paragraph({
    children: runs(line, { size: 24 }),
    spacing: { before: /^\*\*Point/.test(line) ? 320 : 0, after: 200, line: LINE },
    alignment: AlignmentType.LEFT,
  }));
}

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 24, color: "000000" } } } },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children,
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync(OUT, b);
  console.log(`wrote ${OUT} (${(b.length / 1024).toFixed(0)} KB, ${children.length} blocks)`);
});
