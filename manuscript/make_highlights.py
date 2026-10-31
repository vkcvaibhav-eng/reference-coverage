"""Highlights come from the manuscript, so the two can never disagree.

Elsevier allows five bullets of at most 85 characters each.
"""
import io
import sys

LIMIT, MAX_BULLETS = 85, 5

src = io.open("manuscript.md", encoding="utf-8").read()
block = src.split("## Highlights")[1].split("---")[0].strip()
bullets = [b.strip("- ").strip().replace("*", "")
           for b in block.split("\n") if b.strip()]

bad = [b for b in bullets if len(b) > LIMIT]
if bad or len(bullets) > MAX_BULLETS:
    for b in bad:
        print(f"  {len(b)} chars, over {LIMIT}: {b}", file=sys.stderr)
    if len(bullets) > MAX_BULLETS:
        print(f"  {len(bullets)} bullets, over {MAX_BULLETS}", file=sys.stderr)
    sys.exit(1)

io.open("highlights.txt", "w", encoding="utf-8").write("\n".join(bullets) + "\n")
print(f"highlights.txt: {len(bullets)} bullets, "
      f"longest {max(len(b) for b in bullets)} chars")
