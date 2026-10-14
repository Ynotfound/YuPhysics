"""Index user-provided research sources without modifying them."""
import argparse
import hashlib
import json
from pathlib import Path
from pypdf import PdfReader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for path in sorted(args.source.iterdir()):
        if not path.is_file() or path.suffix.lower() not in (".pdf", ".tex"):
            continue
        record = {"path": path.resolve().as_posix(), "name": path.name,
                  "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "bytes": path.stat().st_size}
        if path.suffix.lower() == ".pdf":
            reader = PdfReader(path)
            pages = [page.extract_text() or "" for page in reader.pages]
            text = "\n\n".join(f"--- PDF PAGE {i} ---\n{page}" for i, page in enumerate(pages, 1))
            output = args.output / f"{path.stem}.txt"
            output.write_text(text, encoding="utf-8")
            record.update(pages=len(pages), extracted_text=output.resolve().as_posix())
        records.append(record)
    (args.output / "inventory.json").write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"files": len(records), "pdfs": sum("pages" in r for r in records),
                      "pages": sum(r.get("pages", 0) for r in records),
                      "unique_sha256": len({r["sha256"] for r in records})}))


if __name__ == "__main__":
    main()
