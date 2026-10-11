#!/usr/bin/env python3
"""验证 figure_contract 导出的 PDF/PNG；技术通过不等同人工美观验收。"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import fitz
from PIL import Image, ImageChops


def md5(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "md5").hexdigest()


def validate(manifest_path: Path) -> dict:
    result = {"Manifest": str(manifest_path), "Status": "FAIL", "Errors": ""}
    errors: list[str] = []
    try:
        with manifest_path.open(newline="", encoding="utf-8-sig") as stream:
            records = list(csv.DictReader(stream))
        if len(records) != 1:
            raise ValueError("manifest must contain exactly one figure")
        record = records[0]
        pdf, png = Path(record["PDF"]), Path(record["PNG"])
        result.update(FigureID=record["FigureID"], PDF=str(pdf), PNG=str(png),
                      SourceDataStatus=record["SourceDataStatus"])
        context_fields = ("Project", "Species", "Mode", "AnalysisStage", "AnalysisScope",
                          "Module", "Metric", "Denominator", "Applicability")
        missing_context = [name for name in context_fields if not record.get(name)]
        if missing_context:
            errors.append("missing_context_fields:" + ",".join(missing_context))
        # Explicit unspecified values are honest provenance gaps, not technical
        # rendering failures; retain them visibly in the validation report.
        result.update({name: record.get(name, "") for name in context_fields})
        if md5(pdf) != record["PDFMD5"] or md5(png) != record["PNGMD5"]:
            errors.append("file_hash_changed")
        with fitz.open(pdf) as doc:
            if len(doc) != 1:
                errors.append("expected_one_page_per_figure")
            page = doc[0]
            width_mm, height_mm = page.rect.width * 25.4 / 72, page.rect.height * 25.4 / 72
            result.update(WidthMM=round(width_mm, 3), HeightMM=round(height_mm, 3), PDFPages=len(doc))
            if abs(width_mm - 183) > 0.5 or abs(height_mm - float(record["HeightMM"])) > 0.5:
                errors.append("physical_dimensions")
            fonts = page.get_fonts(full=True)
            result["Fonts"] = ";".join(sorted({font[3] for font in fonts}))
            if any(not doc.extract_font(font[0])[3] for font in fonts):
                errors.append("font_not_embedded")
            if any("LiberationSans" not in font[3].replace(" ", "") for font in fonts):
                errors.append("unexpected_font_family")
            spans = [span for block in page.get_text("dict")["blocks"] if block["type"] == 0
                     for line in block["lines"] for span in line["spans"] if span["text"].strip()]
            if not fonts or not spans:
                errors.append("no_embedded_extractable_text")
            sizes = [span["size"] for span in spans]
            result["MinimumFontPT"] = round(min(sizes), 3) if sizes else ""
            result["MaximumFontPT"] = round(max(sizes), 3) if sizes else ""
            if sizes and (min(sizes) < 5.9 or max(sizes) > 10.2):
                errors.append("font_size_outside_6_10pt")
            outside = sum(s["bbox"][0] < -0.5 or s["bbox"][1] < -0.5 or
                          s["bbox"][2] > page.rect.width + 0.5 or s["bbox"][3] > page.rect.height + 0.5
                          for s in spans)
            result["TextOutsidePage"] = outside
            if outside:
                errors.append("text_outside_page")
            if not spans and not page.get_drawings() and not page.get_images():
                errors.append("empty_pdf")
            expected_px = (round(page.rect.width / 72 * 300), round(page.rect.height / 72 * 300))
        with Image.open(png) as image:
            result.update(PNGWidth=image.width, PNGHeight=image.height)
            if any(abs(a - b) > 1 for a, b in zip(image.size, expected_px)):
                errors.append("png_pixel_dimensions")
            dpi = image.info.get("dpi", (0, 0))
            result["PNGDPI"] = ";".join(str(round(v, 3)) for v in dpi)
            if any(abs(v - 300) > 0.2 for v in dpi):
                errors.append("png_resolution")
            preview = image.convert("RGB")
            preview.thumbnail((512, 512))
            if ImageChops.difference(preview, Image.new("RGB", preview.size, "white")).getbbox() is None:
                errors.append("blank_png")
        source_files = [Path(p) for p in record["SourceDataFiles"].split(";") if p]
        if record["SourceDataStatus"] != "not_supplied" and not source_files:
            errors.append("claimed_source_data_missing")
        if any(not path.is_file() or path.stat().st_size == 0 for path in source_files):
            errors.append("source_data_file_missing")
    except Exception as exc:
        errors.append(f"{type(exc).__name__}: {exc}")
    result["Errors"] = ";".join(errors)
    result["Status"] = "FAIL" if errors else "PASS"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="figure directory, result root, or one .figure_manifest.csv")
    parser.add_argument("--output", type=Path, required=True, help="CSV verification report")
    args = parser.parse_args()
    manifests = [args.root] if args.root.is_file() else sorted(args.root.rglob("*.figure_manifest.csv"))
    if not manifests:
        parser.error("no figure manifests found")
    rows = [validate(path) for path in manifests]
    if args.root.is_dir():
        managed = {Path(row["PDF"]).resolve() for row in rows if row.get("PDF")}
        for pdf in sorted(args.root.rglob("*.pdf")):
            if pdf.resolve() not in managed:
                rows.append({"Manifest": "", "Status": "FAIL", "Errors": "pdf_without_figure_manifest",
                             "FigureID": pdf.stem, "PDF": str(pdf), "PNG": str(pdf.with_suffix(".png"))})
        for png in sorted(args.root.rglob("*.png")):
            if not png.with_suffix(".pdf").is_file():
                rows.append({"Manifest": "", "Status": "FAIL", "Errors": "png_without_pdf_master",
                             "FigureID": png.stem, "PDF": str(png.with_suffix(".pdf")), "PNG": str(png)})
    fields = list(dict.fromkeys(key for row in rows for key in row))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    failures = sum(row["Status"] != "PASS" for row in rows)
    print(f"{len(rows) - failures}/{len(rows)} figures PASS; report: {args.output}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
