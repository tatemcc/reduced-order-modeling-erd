from pathlib import Path

def read_text_file(path):
    return Path(path).read_text(encoding="utf-8")

text = read_text_file("data-relpath.txt")
print(text)
