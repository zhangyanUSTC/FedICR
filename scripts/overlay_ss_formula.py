"""Apply the original formula raster at the original canvas-relative positions."""

import argparse
from PIL import Image


def overlay(chart_path, formula_path, layout, height_fraction=None):
    with Image.open(chart_path) as source:
        chart = source.convert("RGBA")
    with Image.open(formula_path) as source:
        formula = source.convert("RGBA")
    bbox = formula.getchannel("A").getbbox()
    if bbox is None:
        raise ValueError("Formula image has no visible pixels")
    formula = formula.crop(bbox)
    formula = formula.rotate(90, expand=True, resample=Image.Resampling.BICUBIC)
    original_fraction = 0.39 if layout == "single" else 0.222
    fraction = original_fraction if height_fraction is None else height_fraction
    height = round(chart.height * fraction)
    scale = height / formula.height
    formula = formula.resize(
        (max(1, round(formula.width * scale)), height), Image.Resampling.LANCZOS
    )
    x = round(chart.width * 0.681)
    ys = (0.145,) if layout == "single" else (0.088, 0.543)
    for y in ys:
        chart.alpha_composite(formula, (x, round(chart.height * (y + (original_fraction - fraction) / 2))))
    if layout == "combined":
        chart.convert("RGB").save(chart_path, dpi=(300, 300), optimize=True)
    else:
        chart.save(chart_path, dpi=(300, 300))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("chart")
    parser.add_argument("formula")
    parser.add_argument("layout", choices=("single", "combined"))
    parser.add_argument("--height-fraction", type=float)
    args = parser.parse_args()
    overlay(args.chart, args.formula, args.layout, args.height_fraction)
