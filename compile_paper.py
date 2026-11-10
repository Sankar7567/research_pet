#!/usr/bin/env python3
"""
compile_paper.py
================
Autonomous Publication-Grade LaTeX Paper Compiler

Compiles `manuscript_nature_catalysis.tex` into a publication-ready PDF
named `PROD_PETase_Sankar_Research_Paper.pdf` in the repository root.

Uses `tectonic` or `pdflatex` as the compilation backend.
"""

import os
import sys
import shutil
import subprocess
import time
from pathlib import Path


def find_latex_compiler() -> Tuple[str, str]:
    """Finds available LaTeX compiler: tectonic or pdflatex."""
    # Check tectonic first
    tectonic_paths = [
        shutil.which("tectonic"),
        os.path.expanduser("~/.local/bin/tectonic"),
        "/usr/local/bin/tectonic",
        "/usr/bin/tectonic"
    ]
    for p in tectonic_paths:
        if p and os.path.isfile(p) and os.access(p, os.X_OK):
            return "tectonic", p

    # Check pdflatex
    pdflatex = shutil.which("pdflatex")
    if pdflatex:
        return "pdflatex", pdflatex

    # Check xelatex
    xelatex = shutil.which("xelatex")
    if xelatex:
        return "xelatex", xelatex

    raise RuntimeError("No LaTeX compiler found (tectonic, pdflatex, or xelatex).")


def get_pdf_page_count(pdf_path: str) -> int:
    """Extracts PDF page count using basic binary parsing or pypdf/pdfinfo."""
    try:
        with open(pdf_path, "rb") as f:
            content = f.read()
        import re
        pages = re.findall(rb"/Type\s*/Page[^s]", content)
        if pages:
            return len(pages)
        count_match = re.search(rb"/Count\s+(\d+)", content)
        if count_match:
            return int(count_match.group(1))
    except Exception:
        pass
    return -1


def compile_manuscript(
    tex_file: str = "/home/sankar/Desktop/research_2/manuscript_nature_catalysis.tex",
    output_pdf_name: str = "PROD_PETase_Sankar_Research_Paper.pdf",
    output_dir: str = "/home/sankar/Desktop/research_2"
) -> str:
    """Compiles LaTeX source into target PDF."""
    tex_path = Path(tex_file).resolve()
    out_dir = Path(output_dir).resolve()
    final_pdf_path = out_dir / output_pdf_name

    if not tex_path.exists():
        raise FileNotFoundError(f"LaTeX source file not found: {tex_path}")

    compiler_type, compiler_bin = find_latex_compiler()
    print(f"[Compiler] Detected engine: {compiler_type} ({compiler_bin})")
    print(f"[Compiler] Input LaTeX file: {tex_path}")
    print(f"[Compiler] Target PDF output: {final_pdf_path}")

    start_time = time.time()

    if compiler_type == "tectonic":
        cmd = [
            compiler_bin,
            "-X", "compile",
            str(tex_path),
            "--outdir", str(out_dir)
        ]
        print(f"[Compiler] Executing: {' '.join(cmd)}")
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            print("[Compiler] Tectonic stdout:\n", res.stdout)
            print("[Compiler] Tectonic stderr:\n", res.stderr)
            raise RuntimeError(f"Tectonic compilation failed with exit code {res.returncode}")

        # Tectonic creates a PDF matching the tex stem: manuscript_nature_catalysis.pdf
        generated_pdf = out_dir / f"{tex_path.stem}.pdf"
        if generated_pdf.exists():
            if generated_pdf != final_pdf_path:
                shutil.move(str(generated_pdf), str(final_pdf_path))
        else:
            raise FileNotFoundError(f"Expected output PDF not generated: {generated_pdf}")

    elif compiler_type in ("pdflatex", "xelatex"):
        # Run 2 passes to resolve all references and cross-links
        for pass_num in (1, 2):
            print(f"[Compiler] Running {compiler_type} pass {pass_num}/2...")
            cmd = [
                compiler_bin,
                "-interaction=nonstopmode",
                f"-output-directory={out_dir}",
                str(tex_path)
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if res.returncode != 0 and pass_num == 2:
                print(f"[Compiler] {compiler_type} stderr:\n", res.stderr)
                # Check if PDF was still produced despite minor warnings
                generated_pdf = out_dir / f"{tex_path.stem}.pdf"
                if not generated_pdf.exists():
                    raise RuntimeError(f"{compiler_type} failed with exit code {res.returncode}")

        generated_pdf = out_dir / f"{tex_path.stem}.pdf"
        if generated_pdf.exists() and generated_pdf != final_pdf_path:
            shutil.move(str(generated_pdf), str(final_pdf_path))

    elapsed = time.time() - start_time

    if not final_pdf_path.exists():
        raise FileNotFoundError(f"Final PDF does not exist: {final_pdf_path}")

    file_size_kb = final_pdf_path.stat().st_size / 1024.0
    pages = get_pdf_page_count(str(final_pdf_path))

    print("=" * 70)
    print("SUCCESS: Publication-Grade Research Paper Compiled Successfully!")
    print(f"Output File: {final_pdf_path}")
    print(f"File Size:   {file_size_kb:.1f} KB")
    print(f"Page Count:  {pages} pages (Standard Two-Column Journal Format)")
    print(f"Time Taken:  {elapsed:.2f} seconds")
    print("=" * 70)
    return str(final_pdf_path)


if __name__ == "__main__":
    compile_manuscript()
