"""Rewrite the page-number field in a built .docx.

docx-js emits begin/instrText/separate/end inside a single run, which is not
what the format specifies, so this writes the five-run form instead. It also
suppresses line numbering in the footer: the section numbers lines
continuously, the footer contributes one line per page, and the result is a
second number printed in the line-number gutter that tracks the page count.
"""
import re
import shutil
import sys
import zipfile

RPR = ('<w:rPr><w:rFonts w:ascii="Times New Roman" w:cs="Times New Roman" '
       'w:eastAsia="Times New Roman" w:hAnsi="Times New Roman"/>'
       '<w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr>')
FIELD = ('<w:p><w:pPr><w:suppressLineNumbers/><w:jc w:val="center"/></w:pPr>'
         f'<w:r>{RPR}<w:fldChar w:fldCharType="begin"/></w:r>'
         f'<w:r>{RPR}<w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
         f'<w:r>{RPR}<w:fldChar w:fldCharType="separate"/></w:r>'
         f'<w:r>{RPR}<w:t>1</w:t></w:r>'
         f'<w:r>{RPR}<w:fldChar w:fldCharType="end"/></w:r>'
         '</w:p>')

path = sys.argv[1]
tmp = path + ".tmp"
with zipfile.ZipFile(path) as zin, \
        zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
    names = zin.namelist()
    if not any(n.startswith("word/footer") for n in names):
        sys.exit(f"{path}: no footer to fix")
    for name in names:
        data = zin.read(name)
        if name.startswith("word/footer") and name.endswith(".xml"):
            text = data.decode("utf-8")
            text, n = re.subn(r"<w:p>.*?</w:p>", FIELD, text, count=1, flags=re.S)
            if not n:
                sys.exit(f"{name}: no paragraph found")
            data = text.encode("utf-8")
        zout.writestr(name, data)
shutil.move(tmp, path)
print(f"{path}: page-number field rewritten")
