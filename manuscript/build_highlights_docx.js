const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, AlignmentType } = require("docx");
const FONT = "Times New Roman";

const bullets = fs.readFileSync("highlights.txt", "utf8")
  .split(/\r?\n/).map(s => s.trim()).filter(Boolean);

const children = [
  new Paragraph({
    spacing: { after: 240 },
    children: [new TextRun({ text: "Highlights", bold: true, font: FONT, size: 24 })],
  }),
  ...bullets.map(b => new Paragraph({
    bullet: { level: 0 },
    spacing: { after: 120 },
    alignment: AlignmentType.LEFT,
    children: [new TextRun({ text: b, font: FONT, size: 24 })],
  })),
];

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 24 } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 },
    margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } }, children }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync("highlights.docx", b);
  console.log(`wrote highlights.docx (${bullets.length} bullets)`);
});
